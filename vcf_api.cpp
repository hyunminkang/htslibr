#include<Rcpp.h>
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
using namespace Rcpp;
using namespace std;

//' extract values from the INFO field
//' @param vcf the VCF/BCF file path
//' @param index the CSI/TBI index file path
//' @param reg a region query of the form: chr:start-end 
//' @param tag the field in the INFO field to extract. Only accepts one string value at this time. Only can extract numeric fields at the moment. 
//' @description Use this function to extract the INFO field values for a single INFO field in a give
//' region based query. 
//' @return a dataframe with the chrom, pos, and value of the given INFO field
//' @examples
//' \dontrun{extract_info(vcf, index, "1:10001-100500", "AC")}
// [[Rcpp::export]]
DataFrame extract_info(std::string vcf, std::string index, std::string& reg, std::string &tag) {
    htsFile *fp = hts_open(vcf.c_str(), "r");
    if (!fp) REprintf("couldn't read vcf %s", vcf);
    const htsFormat* fmt = &fp->format;
    char* description = hts_format_description(fmt);
    std::string format(description);
    Rcout << "detecting format " << format << endl;
    bool use_csi = false;
    if (index.find("csi") != std::string::npos) {
        use_csi = true;
        Rcout << "using csi index" << endl;
    }
    hts_idx_t *csi_idx;
    tbx_t *tbi_idx;
    if (use_csi) {
         csi_idx = bcf_index_load2(vcf.c_str(), index.c_str());
        if(!csi_idx) {
            stop("csi index is null");
        }
    } else {
        tbi_idx = tbx_index_load2(vcf.c_str(), index.c_str());
        if(!tbi_idx) {
            stop("tbi index is null");
        }
    }

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) REprintf("can't reader header for vcf %s", vcf);
    bcf1_t *line = bcf_init();

    hts_itr_t *itr;
    if (use_csi) {
        itr = bcf_itr_querys(csi_idx, hdr, reg.c_str());
    } else {
        itr = tbx_itr_querys(tbi_idx, reg.c_str());
        if(!itr) stop("itr is null");
    }
    int r;
    int ret;
    kstring_t s = {0, 0, NULL}; // needs to be initialized to prevent segfault    
    bcf_info_t *info;

    // report chrom and position
    CharacterVector chroms;
    IntegerVector positions;

    // don't know what type the INFO field 'tag' will have so we initialize three possible types
    IntegerVector int_res;
    NumericVector num_res;
    CharacterVector char_res;
    bool is_int = false;
    bool is_float = false;
    bool is_string = false;

    if (use_csi) {
        // TODO: refactor so there is less repeated code
        while((r = bcf_itr_next(fp, itr, line)) >= 0) {

            chroms.push_back(bcf_hdr_id2name(hdr, line->rid));
            positions.push_back(line->pos);

            info = bcf_get_info(hdr, line, tag.c_str());
            if (info == nullptr) stop("info field does not exist");

            if (info->len == 1) {
                    // see here as a reference https://github.com/brentp/cyvcf2/blob/master/cyvcf2/cyvcf2.pyx#L2003
                if (info->type == BCF_BT_INT8 || info->type == BCF_BT_INT16 || info->type == BCF_BT_INT32) {
                    int_res.push_back(info->v1.i); 
                    is_int = true;
                } else if(info->type == BCF_BT_FLOAT) {
                    // num_res.push_back(bcf_float_is_missing(info->v1.f)); // appears to convert values to 0
                    num_res.push_back(info->v1.f);
                    is_float = true;
                }
            }
        }
    } else {
        while((r = tbx_itr_next(fp, tbi_idx, itr, &s)) >= 0) {

            ret = vcf_parse(&s, hdr, line);
            if(ret > 0) stop("vcf parsing error");

            chroms.push_back(bcf_hdr_id2name(hdr, line->rid));
            positions.push_back(line->pos);

            info = bcf_get_info(hdr, line, tag.c_str());
            if (info == nullptr) stop("info field does not exist");

            if (info->len == 1) {
                if (info->type == BCF_BT_INT8 || info->type == BCF_BT_INT16 || info->type == BCF_BT_INT32) {
                    int_res.push_back(info->v1.i); 
                    is_int = true; // need this for the final return type
                } else if(info->type == BCF_BT_FLOAT) {
                    // num_res.push_back(bcf_float_is_missing(info->v1.f)); // appears to convert values to 0
                    num_res.push_back(info->v1.f);
                    is_float = true;
                }
            }
        }
    }

    free(s.s);

    bcf_destroy(line);
    hts_close(fp);

    if (is_int) {
        return DataFrame::create(
            Named("chrom") = chroms,
            Named("pos") = positions,
            Named("value") = int_res
        );
    } else if(is_float) {
        return DataFrame::create(
            Named("chrom") = chroms,
            Named("pos") = positions,
            Named("value") = num_res
        );
    } else {
        return DataFrame::create(
            Named("chrom") = chroms,
            Named("pos") = positions,
            Named("value") = char_res
        );
    }

    stop("not yet implemented"); // this happens if the INFO field is a string
}

//' extract the genotypes for a given region from the GT field
//' @param vcf the VCF/BCF file path
//' @param index the CSI/TBI index file path
//' @param reg a region query of the form: chr:start-end 
//' @description Use this function to extract the genotypes from the GT field. Will return as a 
//' IntegerMatrix of dimensions haplotypes x variants. That is, each (diploid) individual will have two consecutve rows.
//' No existing support for using the phase of the genotypes (if present) or for handling missing values or
//' variable ploidy. 
//' @return a integer matrix of dimension (number of haplotypes x number of variants).
//' @examples
//' \dontrun{extract_genotypes(vcf, index, "1:10001-100500")}
// [[Rcpp::export]]
SEXP extract_genotypes(std::string vcf, std::string index, std::string& reg) {
    htsFile *fp = hts_open(vcf.c_str(), "r");
    if (!fp) REprintf("couldn't read vcf %s", vcf);
    const htsFormat* fmt = &fp->format;
    char* description = hts_format_description(fmt);
    std::string format(description);
    Rcout << "detecting format " << format << endl;
    bool use_csi = false;
    if (index.find("csi") != std::string::npos) {
        use_csi = true;
        Rcout << "using csi index" << endl;
    }
    hts_idx_t *csi_idx;
    tbx_t *tbi_idx;
    if (use_csi) {
         csi_idx = bcf_index_load2(vcf.c_str(), index.c_str());
        if(!csi_idx) {
            stop("csi index is null");
        }
    } else {
        tbi_idx = tbx_index_load2(vcf.c_str(), index.c_str());
        if(!tbi_idx) {
            stop("tbi index is null");
        }
    }

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) REprintf("can't reader header for vcf %s", vcf);
    bcf1_t *line = bcf_init();

    hts_itr_t *itr;
    if (use_csi) {
        itr = bcf_itr_querys(csi_idx, hdr, reg.c_str());
    } else {
        itr = tbx_itr_querys(tbi_idx, reg.c_str());
        if(!itr) stop("itr is null");
    }
    int r;
    int ret;
    kstring_t s = {0, 0, NULL}; // needs to be initialized to prevent segfault    

    bcf_info_t *info;
    int i, j, ngt, n_samples = bcf_hdr_nsamples(hdr);
    int32_t *gt_arr = NULL, ngt_arr = 0;
    int num_variants = 0;
    // int32_t **gt_all = NULL;
    // std::vector<IntegerVector> genotype_matrix; // rows are variants, columns are genotypes
    IntegerVector genotypes;

    Rprintf("detecting %d samples\n", n_samples);
    if (use_csi) {
        // TODO: refactor so there is less repeated code
        while((r = bcf_itr_next(fp, itr, line)) >= 0) {
            num_variants++;

            // https://github.com/samtools/htslib/blob/2da4c7dd951428fa9d0d4049394045d1cace4133/htslib/vcf.h#L790
            ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
            int max_ploidy = ngt / n_samples;
            if(max_ploidy != 2) stop("currently only support for diploid organisms");

            for (i = 0; i < n_samples; i++) {
                int32_t *ptr = gt_arr + i * max_ploidy; // pointer to ith individual's genotypes

                for(j = 0; j < max_ploidy; j++) {

                    if ( ptr[j]==bcf_int32_vector_end ) break;

                    if ( bcf_gt_is_missing(ptr[j]) ) stop("missing value support not yet added");

                    int allele_index = bcf_gt_allele(ptr[j]);
                    // Rprintf("sample i: %d and allele_index %d\n", i, allele_index);
                    genotypes.push_back(allele_index);
                }
            }

        }
    } else {
        while((r = tbx_itr_next(fp, tbi_idx, itr, &s)) >= 0) {

            ret = vcf_parse(&s, hdr, line);
            if(ret > 0) stop("vcf parsing error");

            num_variants++;

            ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
            int max_ploidy = ngt / n_samples;
            if(max_ploidy != 2) stop("currently only support for diploid organisms");

            for (i = 0; i < n_samples; i++) {
                int32_t *ptr = gt_arr + i * max_ploidy; // pointer to ith individual's genotypes

                for(j = 0; j < max_ploidy; j++) {

                    if ( ptr[j]==bcf_int32_vector_end ) break;

                    if ( bcf_gt_is_missing(ptr[j]) ) stop("missing value support not yet added");

                    int allele_index = bcf_gt_allele(ptr[j]);
                    genotypes.push_back(allele_index);
                }
            }

        }
    }


    bcf_destroy(line);
    hts_close(fp);

    // got this idea from https://stackoverflow.com/questions/19864226/convert-stdvector-to-rcpp-matrix
    genotypes.attr("dim") = Dimension(ngt, num_variants); // haplotypes x variants

    return genotypes;
}
