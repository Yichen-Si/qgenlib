#ifndef __SEQ_UTILS_H
#define __SEQ_UTILS_H

#include <cstdint>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <functional>
#include <vector>
#include "htslib/sam.h"
#include "htslib/vcf.h"

#define MAX_NT5_UNIT_64 27
#define MAX_NT5_UNIT_128 55

typedef __uint128_t uint128_t; // 128-bit unsigned integer
// custom hash function for uint128_t
struct uint128_hash {
    std::size_t operator()(const uint128_t& key) const {
        return std::hash<uint64_t>()(key >> 64) ^ std::hash<uint64_t>()(key);
    }
};

extern const unsigned char seq_nt6_table[256];
extern const unsigned char seq_nt5_table[256];
extern const char *seq_nt16_rev_table;
extern const char *seq_nt5_rev_table;
extern const char *seq_nt4_rev_table;
extern const unsigned char seq_nt16to4_table[16];
extern const unsigned char seq_nt16comp_table[16];
extern const int bitcnt_table[16];
extern const char comp_tab[128];

int32_t seq_iupac_mismatch(const char* seq, const char* pattern, int32_t len);
void seq_revcomp(char* seq, int32_t len);
uint64_t seq2bits(const char* seq, int32_t len, uint8_t nonACGTs = 0x03);
uint64_t seq2nt5(const char* seq, int32_t len);
bool nt52seq(uint64_t nt5, int32_t len, char *seq);
int32_t seq2nt5multi(const char *seq, int32_t lseq, uint64_t *nt5s, int32_t nt5unit = 27);
bool nt5multi2seq(uint64_t *nt5s, int32_t lseq, char *seq, int32_t nt5unit = 27);
uint128_t seq2bits_128(const char *seq, int32_t len, uint8_t nonACGTs = 0x03);
uint128_t seq2nt5_128(const char *seq, int32_t len);
bool nt52seq_128(uint128_t nt5, int32_t len, char *seq);

int32_t nt4_hamming_dist(uint64_t a, uint64_t b, int32_t max_dist = 32);
uint64_t seq2bits2(const char* seq, int32_t len, std::vector<uint8_t>& nonACGTs, uint8_t ambig_force = 1);

// utilities for SAM/BAM/CRAM files
int32_t bam_flex_read(samFile* in, hts_itr_t* iter, bam_hdr_t* hdr, bam1_t* b);

// utilities for VCF/BCF files
int32_t bcf_flex_read(vcfFile* in, hts_itr_t* iter, bcf_hdr_t* hdr, bcf1_t* b);
//const char* bcf_get_chrom(bcf_hdr_t *h, bcf1_t *v);
bool extract_ac_an_af_from_info(bcf1_t* b, bcf_hdr_t* hdr, int32_t* p_ac, int32_t* p_an, double* p_af);

#endif // __SEQ_UTILS_H
