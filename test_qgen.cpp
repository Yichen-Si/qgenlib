#include "qgenlib/params.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_utils.h"
#include "qgenlib/hts_utils.h"

#include <cstdio>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include <set>
#include <ctime>
#include <cmath>
#include <cassert>

#include <map>

struct AlleleStats {
    int count = 0;
    int total_bq = 0;
    int min_bq = 255;
    int max_bq = 0;
    int forward = 0;
    int backward = 0;
    std::vector<int> sb_list;
    std::vector<int> bq_list;
};

int32_t main(int32_t argc, char** argv) {
  std::string inPileup, output, outpref;
  int32_t min_ad = 1, min_bq = 0, max_bq = 30, mean_bq = 0;
  int32_t c_base = 4, c_bq = 5, c_sb = -1;
  int32_t debug = 0, verbose = 1000;
  std::vector<std::string> colNames;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input", NULL)
    LONG_STRING_PARAM("pileup",&inPileup, "Input pileup")
    LONG_MULTI_STRING_PARAM("columns",&colNames, "Columns after the default BQ")
    LONG_INT_PARAM("col-sb",&c_sb, "If present, which column contains spatial barcodes")
    LONG_INT_PARAM("debug",&debug, "Debug")

    LONG_PARAM_GROUP("Options for output", NULL)
    LONG_STRING_PARAM("outpref",&outpref, "Output prefix")
    LONG_INT_PARAM("min-ad",&min_ad, "Minimum AD for ALT")
    LONG_INT_PARAM("min-bq",&min_bq, "Minimum min BQ for ALT")
    LONG_INT_PARAM("max-bq",&max_bq, "Minimum max BQ for ALT")
    LONG_INT_PARAM("mean-bq",&mean_bq, "Minimum mean BQ for ALT")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( inPileup.empty() || outpref.empty() )
    error("[E:%s:%d %s] --input --output are required but at least one is missing",__FILE__,__LINE__,__PRETTY_FUNCTION__);

  if (c_sb < 0) {
    for (uint32_t i=0; i < colNames.size(); ++i) {
      if ( colNames[i] == "SB" ) c_sb = i;
    }
  }
  int32_t n_cols = 6;
  bool output_per_sb = (c_sb >= 0);
  if (output_per_sb)
    n_cols = std::max(n_cols, c_sb + 1);

  output = outpref + ".snv.summary.tsv.gz";
  htsFile* wfs = hts_open(output.c_str(), "wz");
  hprintf(wfs, "#CHROM\tPOS\tREF\tALLELE\tDP\tAD\tFwd\tRev\tminBQ\tmaxBQ\tmeanBQ\n");
  htsFile* wf = NULL;

  if (wfs == NULL)
    error("[E:%s:%d %s] Cannot open %s for writing",__FILE__,__LINE__,__PRETTY_FUNCTION__,output.c_str());

  if (output_per_sb) {
    output = outpref + ".snv.tsv.gz";
    wf = hts_open(output.c_str(), "wz");
    hprintf(wf, "#CHROM\tPOS\tREF\tALT\tSTRAND\tBQ\tSB\n");
  }

  tsv_reader tr(inPileup.c_str());
  int32_t n_skip = 0, n_snv = 0;
  while(tr.read_line() > 0) {
    if (tr.nfields < n_cols) {
      error("[E:%s:%d %s] At the %llu line the input file does not have enough fields",__FILE__,__LINE__,__PRETTY_FUNCTION__,tr.nlines);
    }
    if (tr.nlines % verbose == 0) {
      notice("Processed %llu lines, skipped %d, output %d SNVs", tr.nlines, n_skip, n_snv);
    }

    int32_t nreads = tr.int_field_at(3);
    char ref_allele = *tr.str_field_at(2);
    char ref_reverse = std::tolower(ref_allele);
    std::string bases(tr.str_field_at(c_base));
    const char* bqs = tr.str_field_at(c_bq);
    const char* sbs = output_per_sb ? tr.str_field_at(c_sb) : NULL;
    std::vector<std::string> sb_list;
    if (output_per_sb) {
      split(sb_list, ",", sbs);
    }
    // find SNV, ignoring indels for now
    std::map<char, AlleleStats> allele_map;
    int pt_read = 0;
    std::map<char, int32_t> base_ct; // For debug
    for (size_t i = 0, j = 0; i < bases.size();  ++i) {
        // char base = std::toupper(bases[i]);
        char base = bases[i];
        char base_abs = std::toupper(base);
        base_ct[base_abs]++;
// if (debug && (base == '*' || base == '#' || base == '>' || base == '<')) {
//   std::cout << base << '\t' << bqs[j] << '\t' << i << '/' << bases.size() << '\t' << j << '/' << nreads << '\n';
// }
        if (base == '.') {
          base = ref_allele;
        } else if (base == ',') {
          base = ref_reverse;
        } else if (base == '^') { // First base in the read, followed by MQ
          i++; // Skip the next character (MQ)
          continue;
        } else if (base == '$') { // Last base in the read
          continue;
        } else if (base == '+' || base == '-') { // Indels
          size_t pt_end = bases.find_first_not_of("0123456789", i + 1);
          int indel_len = std::stoi(bases.substr(i + 1, pt_end - i - 1));
          i = pt_end + indel_len - 1; // The last base of the indel
          // if (base == '+')
          //   j += indel_len; // Skip the insertion in the base quality string
          continue;
        } else if (base == '*' || base == '#') { // Deletion
          pt_read++;
          j++; // Mapped to indel?
          continue;
        } else if (base == '>' || base == '<') { // Reference skip
          pt_read++;
          j++;
          continue;
        }
        // It is a SNV or REF
        int quality = bqs[j] - 33; // Convert to Phred scale
if (quality < 0) { // shouldn't happen
  std::cout << "Quality is negative: [" << tr.int_field_at(1) << ']' << bqs[j] << ',' << base << ',' << quality << ',' << i << '/' << bases.size() << ", " << j << '(' << pt_read << '/' << nreads << ')' << std::endl;
  quality = 0;
  for (const auto& kv : base_ct) {
    std::cout << kv.first << ':' << kv.second << '\t';
  }
  std::cout << std::endl;
}
        base_abs = std::toupper(base);
        allele_map[base_abs].count++;
        allele_map[base_abs].total_bq += quality;
        allele_map[base_abs].min_bq = std::min(allele_map[base_abs].min_bq, quality);
        allele_map[base_abs].max_bq = std::max(allele_map[base_abs].max_bq, quality);
        if (output_per_sb) {
          if (base < 97) {
            allele_map[base_abs].sb_list.push_back(pt_read+1);
          } else {
            allele_map[base_abs].sb_list.push_back(-pt_read-1);
          }
          allele_map[base_abs].bq_list.push_back(quality);
        }
        if (base < 97) {
          allele_map[base_abs].forward++;
        } else {
          allele_map[base_abs].backward++;
        }

        pt_read++;
        j++;

if (debug && i == bases.size() - 1) {
  // if (nreads != pt_read || bqs[j] != '\0') {
    std::string bqstr(bqs);
    std::cout << '[' << tr.int_field_at(1) << ']' << nreads << '\t' << i << '\t' << j << '\t' << bqs[j] << '\t' << bqstr.length() << '\t' << pt_read << '\t' << sb_list.size() << std::endl;
    for (const auto& kv : base_ct) {
      std::cout << kv.first << ':' << kv.second << '\t';
    }
    std::cout << std::endl;
  // }
}

    }

    std::vector<char> kept_allele;
    if (allele_map.size() > 1) {
      for (const auto& kv : allele_map) {
        int32_t mean_bq_rd = (int32_t) std::round(kv.second.total_bq / kv.second.count);

if (debug && debug >= n_snv) {
  std::cout << tr.int_field_at(1) << "\t" << ref_allele
  << "\t" << kv.first << "\t"
  << kv.second.count << '(' << kv.second.forward << ',' << kv.second.backward << ')'
  << "\t" << kv.second.min_bq << "\t" << kv.second.max_bq << "\t" << mean_bq_rd << "\t";
  if (output_per_sb) {
    std::cout << kv.second.sb_list.size() << "," << kv.second.bq_list.size() << std::endl;
  } else {
    std::cout << std::endl;
  }
}

        if (kv.second.count < min_ad) continue;
        if (kv.second.min_bq < min_bq) continue;
        if (kv.second.max_bq < max_bq) continue;

        if (mean_bq_rd < mean_bq) continue;
        kept_allele.push_back(kv.first);
      }
    }
    if (kept_allele.size() < 2) {
      n_skip++;
      continue;
    }

    // Output
    for (const auto & k : kept_allele) {
      AlleleStats& allele = allele_map[k];
      int32_t mean_bq_rd = (int32_t) std::round(allele.total_bq / allele.count);

      // "#CHROM\tPOS\tREF\tALLELE\tDP\tAD\tFwd\tRev\tminBQ\tmaxBQ\tmeanBQ\n"
      hprintf(wfs, "%s\t%d\t%c\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", tr.str_field_at(0), tr.int_field_at(1), ref_allele, k, nreads, allele.count, allele.forward, allele.backward, allele.min_bq, allele.max_bq, mean_bq_rd);

      if (output_per_sb) {
        // "#CHROM\tPOS\tREF\tALT\tSTRAND\tBQ\tSB\n"
        for (size_t i = 0; i < allele.sb_list.size(); ++i) {
          int sb_id = allele.sb_list[i];
          char strand = '+';
          if (sb_id < 0) {
            sb_id = - sb_id - 1;
            strand = '-';
          } else {
            sb_id = sb_id - 1;
          }
          hprintf(wf, "%s\t%d\t%c\t%c\t%c\t%d\t%s\n", tr.str_field_at(0), tr.int_field_at(1), ref_allele, k, strand, allele.bq_list[i], sb_list[sb_id].c_str() );
        }
      }
    }
    n_snv++;

if (debug && debug < n_snv) {break;}

  }

  hts_close(wfs);
  if (wf != NULL) hts_close(wf);

  return 0;
}
