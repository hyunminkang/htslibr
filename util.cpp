#include<Rcpp.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
using namespace Rcpp;
using namespace std;

//' title print htslib version
//' Print out htslib version
// [[Rcpp::export]]
void htslib_version() {
    Rcout << hts_version() << endl;
}

//' Detect format of the file
//' @param fname the file to detect the format of
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
//[[Rcpp::export]]
CharacterVector extract_sequence(std::string bam, std::string index, std::string reg) {
    htsFile *fp = hts_open(bam.c_str(), "r");
    hts_idx_t *idx = sam_index_load(fp, index.c_str());
    // bam_hdr_t *hdr = NULL;
    // hdr = sam_hdr_read(fp);
    bam_hdr_t *hdr = sam_hdr_read(fp);
    // Rcout << hdr->text << endl;

    // int beg;
    // int end;
    // Rcout << "1" << endl;
    // const char* reg = hts_parse_reg("chr1:10001-10026", &beg, &end);
    // Rcout << "reg is" << reg << endl;
    
    // Rcout << hdr << endl;
    // hts_itr_t *itr = sam_itr_querys(idx, hdr, "chr1:10001-10026");
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
        Rcout << "target id: " << itr->tid << endl;
        Rcout << "target name: " << hdr->target_name[itr->tid] << endl;
        c = &b->core;
        seq = bam_get_seq(b);
        seq_str = "";
        for(int i = 0; i < c->l_qseq; ++i) {
            seq_str += seq_nt16_str[bam_seqi(seq, i)];
        }
        // Rcout << "seq is: " << seq_str << endl;
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

//' Estimate the depth for each position in a given region
//' @param bam the cram/bam/sam file
//' @param index the index of the cram/bam/sam file
//' @param reg the region of interest, typically in format of chr1:start-begin
//[[Rcpp::export]]
DataFrame depth(std::string bam, std::string index, const std::string& reg, const int flank_bp) {
    htsFile *fp = hts_open(bam.c_str(), "r");

    hts_idx_t *idx = sam_index_load(fp, index.c_str());
    bam_hdr_t *hdr = sam_hdr_read(fp);

    int beg_reg;
    int end_reg;
    const char* reg_parsed = hts_parse_reg(reg.c_str(), &beg_reg, &end_reg);
    std::string tid_name("");
    const char* separator = ":";

    Rcout << "reg parsed " << reg_parsed << endl;

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
    // + 1 because hts_parse_reg appears to subract 1 from beg
    int beg_reg_flank = beg_reg - flank_bp  + 1 > 0 ? beg_reg - flank_bp  + 1: 0; 
    int end_reg_flank = end_reg + flank_bp < target_len ? end_reg + flank_bp : target_len;

    Rcout << "beg_reg " << beg_reg << " and beg_reg_flank" << beg_reg_flank << endl;

    std::string flank_region = tid_name + std::string(":") + std::to_string(beg_reg_flank) + std::string("-") + std::to_string(end_reg_flank);
    
    Rcout << "flank region " << flank_region << endl;
    hts_itr_t *itr = sam_itr_querys(idx, hdr, flank_region.c_str());
    hts_idx_destroy(idx);

    IntegerVector coverage(target_len);
    CharacterVector chrom(target_len);

    bam1_t *b = NULL;
    b = bam_init1();
    int r = 0;
    bam1_core_t *c = NULL;
    int start = 0;
    int end = 0;
    while((r = sam_itr_next(fp, itr, b)) >= 0) {
        c = &b->core;
        start =  c->pos;
        end =  bam_endpos(b);
        Rcout << "start " << start << " end " << end << endl;
        coverage[start]++; // incremement start of read
        coverage[end]--; // decrement end of read
    }
    hts_itr_destroy(itr);

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
