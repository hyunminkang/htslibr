#include<Rcpp.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
using namespace Rcpp;
using namespace std;

//' @title print htslib version
//' @description Print out htslib version
// [[Rcpp::export]]
void htslib_version() {
    Rcout << hts_version() << endl;
}

//' Detect format of the file
//' @param fname the file to detect the format of
//' @return a string with the file format description
// [[Rcpp::export]]
std::string check_format(std::string fname) {
    htsFile *fp = hts_open(fname.c_str(), "r");
    const htsFormat *fmt = hts_get_format(fp);
    return std::string(hts_format_description(fmt));
}

//' Extract the sequences for a given region
//' @param bam the cram/bam/sam file
//' @param index the index of the cram/bam/sam file
//' @param reg the region of interest, typically in format of chr1:start-begin
//' @return a character vector with the sequences in the given region
//' @examples
//' \dontrun{count_kmer(bam, index, "chr1:10001-100050")}
//[[Rcpp::export]]
CharacterVector extract_sequence(std::string bam, std::string index, std::string reg) {
    htsFile *fp = hts_open(bam.c_str(), "r");
    hts_idx_t *idx = sam_index_load(fp, index.c_str());
    bam_hdr_t *hdr = sam_hdr_read(fp);

    
    hts_itr_t *itr = sam_itr_querys(idx, hdr, reg.c_str());
    hts_idx_destroy(idx);

    bam1_t *b = NULL;
    b = bam_init1();
    uint8_t *seq = NULL;
    int r = 0;
    bam1_core_t *c = NULL;
    std::string seq_str("");
    CharacterVector sequences;
    while((r = sam_itr_next(fp, itr, b)) >= 0) {
        c = &b->core;
        seq = bam_get_seq(b);
        seq_str = "";
        for(int i = 0; i < c->l_qseq; ++i) {
            seq_str += seq_nt16_str[bam_seqi(seq, i)];
        }
        sequences.push_back(seq_str);
    }
    hts_itr_destroy(itr);
    return sequences;
}

// reference: https://rosettacode.org/wiki/Count_occurrences_of_a_tubstring#C.2B.2B
int count_kmer_seq(const std::string& seq, const std::string& kmer)
{
    if (kmer.length() == 0) return 0;
    int count = 0;
    for (size_t offset = seq.find(kmer); offset != std::string::npos;
            offset = seq.find(kmer, offset + kmer.length()))
    {
        ++count;
    }
    return count;
}

//' count the number of times a kmer is present in a region
//' @param bam the cram/bam/sam file
//' @param index the index of the cram/bam/sam file
//' @param reg the region of interest, typically in format of chr1:start-begin
//' @param kmer the substring to search for in the reads
//' @return a dataframe with the sequnce reads and counts of the given kmer per read (i.e. two columns)
//' @examples
//' \dontrun{count_kmer(bam, index, "chr1:10001-100050", "TTACGG")}
// [[Rcpp::export]]
DataFrame count_kmer(std::string bam, std::string index, const std::string& reg, const std::string& kmer) {
    htsFile *fp = hts_open(bam.c_str(), "r");

    int count = 0;
    IntegerVector counts; 
    hts_idx_t *idx = sam_index_load(fp, index.c_str());
    bam_hdr_t *hdr = sam_hdr_read(fp);
    hts_itr_t *itr = sam_itr_querys(idx, hdr, reg.c_str());
    hts_idx_destroy(idx);

    bam1_t *b = NULL;
    b = bam_init1();
    int r = 0;
    uint8_t *seq = NULL;
    bam1_core_t *c = NULL;
    std::string seq_str("");
    CharacterVector sequences;
    while((r = sam_itr_next(fp, itr, b)) >= 0) {
        c = &b->core;
        seq = bam_get_seq(b);
        seq_str = "";
        for(int i = 0; i < c->l_qseq; ++i) {
            seq_str += seq_nt16_str[bam_seqi(seq, i)];
        }
        sequences.push_back(seq_str);
        count = count_kmer_seq(seq_str, kmer);
        counts.push_back(count);
    }
    hts_itr_destroy(itr);

    return DataFrame::create(
        Named("seq") = sequences,
        Named("counts") = counts
    );
}

//' Calculate the GC content for a region
//' @param bam the cram/bam/sam file
//' @param index the index of the cram/bam/sam file
//' @param reg the region of interest, typically in format of chr1:start-begin
//' @return a dataframe with the sequnce reads, counts of GC bases, and proportion of GC per read
//' @examples
//' \dontrun{gc_content(bam, index, "chr1:10001-100050")}
//[[Rcpp::export]]
DataFrame gc_content(std::string bam, std::string index, const std::string& reg) {
    htsFile *fp = hts_open(bam.c_str(), "r");

    IntegerVector counts; 
    NumericVector props; 
    hts_idx_t *idx = sam_index_load(fp, index.c_str());
    bam_hdr_t *hdr = sam_hdr_read(fp);
    hts_itr_t *itr = sam_itr_querys(idx, hdr, reg.c_str());
    hts_idx_destroy(idx);

    bam1_t *b = NULL;
    b = bam_init1();
    uint8_t *seq = NULL;
    int r = 0;
    bam1_core_t *c = NULL;
    std::string seq_str("");
    CharacterVector sequences;
    int count_c = 0;
    int count_g = 0;
    int count_gc = 0;
    while((r = sam_itr_next(fp, itr, b)) >= 0) {
        c = &b->core;
        seq = bam_get_seq(b);
        seq_str = "";
        for(int i = 0; i < c->l_qseq; ++i) {
            seq_str += seq_nt16_str[bam_seqi(seq, i)];
        }
        sequences.push_back(seq_str);
        count_c = count_kmer_seq(seq_str, "C");
        count_g = count_kmer_seq(seq_str, "G");
        count_gc = count_c + count_g;
        counts.push_back(count_gc);
        props.push_back(count_gc * 1.0 / seq_str.length());
    }
    hts_itr_destroy(itr);

    return DataFrame::create(
        Named("seq") = sequences,
        Named("gc_count") = counts,
        Named("gc_prop") = props
    );
}

//' Estimate approximate depth for each position in a given region
//' @param bam the cram/bam/sam file
//' @param index the index of the cram/bam/sam file
//' @param reg the region of interest, typically in format of chr1:start-begin
//' @description Calculate an approximate measure for a given region in a CRAM/BAM file. 
//' @details This is only an approximate depth, based on the 'fast mode' algorithm from mosdepth.
//' It allocates an array for the entire chromosome, an increments each start site and decrements
//' each end site, and then takes the cumulative sum. It does not account for mismatches in the alignment,
//' hence, the reference to an 'approximate' depth. 
//' @return a dataframe with the chrom, position, and approximate depth (three columns)
//' @examples
//' \dontrun{depth(bam, index, "chr1:10001-100050")}
//[[Rcpp::export]]
DataFrame depth(std::string bam, std::string index, const std::string& reg) {
    htsFile *fp = hts_open(bam.c_str(), "r");

    hts_idx_t *idx = sam_index_load(fp, index.c_str());
    bam_hdr_t *hdr = sam_hdr_read(fp);

    int beg_reg;
    int end_reg;
    const char* reg_parsed = hts_parse_reg(reg.c_str(), &beg_reg, &end_reg);
    std::string tid_name("");
    const char* separator = ":";

    // parse out the chrom from the region
    for(int i = 0; i < reg.find(separator); i++) {
        tid_name += reg[i]; 
    }

    int n_targets = hdr->n_targets;
    char **target_names = hdr->target_name;

    int tid = 0;
    for(tid = 0; tid < n_targets; tid++) {
        std::string target_name = target_names[tid]; // find tid corresponding to chrom in region
        if(target_name.compare(tid_name) == 0) {
            break;
        }
    }

    // get length of chromosome
    int target_len = hdr->target_len[tid];
    
    hts_itr_t *itr = sam_itr_querys(idx, hdr, reg.c_str());

    IntegerVector coverage(target_len);
    CharacterVector chrom(target_len);

    bam1_t *b = NULL;
    b = bam_init1();
    int r = 0;
    bam1_core_t *c = NULL;
    int start = 0;
    int end = 0;
    while((r = sam_itr_next(fp, itr, b)) >= 0) {
        // see https://github.com/brentp/mosdepth/blob/master/mosdepth.nim#L308
        // similar to the mosdepth fast algorithm
        c = &b->core;
        start =  c->pos;
        end =  bam_endpos(b);
        coverage[start]++; // incremement start of read
        coverage[end]--; // decrement end of read
    }
    hts_itr_destroy(itr);
    hts_idx_destroy(idx);

    IntegerVector depths = cumsum(coverage);

    IntegerVector pos(target_len);
    
    std::iota(pos.begin(), pos.end(), 0);
    std::fill(chrom.begin(), chrom.end(), tid_name);

    IntegerVector reg_indices(end_reg - beg_reg);

    std::iota(reg_indices.begin(), reg_indices.end(), beg_reg);

    return DataFrame::create(
        Named("chrom") = chrom[reg_indices],
        Named("pos") = pos[reg_indices],
        Named("depth") = depths[reg_indices]
    );
}

