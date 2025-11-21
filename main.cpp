#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <zlib.h>

#include "cnpy.h"
#include "kmc_api/kmc_file.h"
#include "kseq.h"
extern "C" {
#include "kmer.h"
}

KSEQ_INIT(gzFile, gzread)

static const char *const USAGE_MESSAGE =
    "Usage: fcgr [-m MASK] <in-kmc-list>\n"
    "Options:\n"
    "        -m MASK       use this mask (default: 1{k})\n"
    "        -o DIR        store output(s) in this directory (default: ./)\n"
    "        -p            add mask to output filenames (default: false)\n"
    "        -h            display this help and exit\n";

/**
// Just for debugging
std::string get_skmer(const uint64_t &kmer, const uint32_t &klen) {
  std::string kmer_s(klen, 'X');
  for (uint32_t i = 0; i < klen; ++i) {
    kmer_s[klen - i - 1] = "ACGT"[(kmer >> (2 * i)) & 3];
  }
  return kmer_s;
}
**/

void fill_index(std::vector<uint64_t> &index, uint32_t k) {
  double nt_coordX[] = {
      1.0,
      -1.0,
      -1.0,
      1.0,
  };
  double nt_coordY[] = {
      1.0,
      1.0,
      -1.0,
      -1.0,
  };
  uint32_t hl = (1 << k);
  uint64_t kmer = 0;
  double x = 0.0, y = 0.0;
  uint64_t i = 0, j = 0, p = 0;
  for (kmer = 0; kmer < (1UL << k * 2); ++kmer) {
    x = y = 0;
    for (int _ = k - 1; _ >= 0; --_) {
      x = (x + nt_coordX[(kmer >> (_ * 2)) & 3]) / 2.0;
      y = (y + nt_coordY[(kmer >> (_ * 2)) & 3]) / 2.0;
    }
    x = (x + 1) / 2;
    y = (y + 1) / 2;

    i = (hl - ceil(y * hl) + 1) - 1;
    j = ceil(x * hl) - 1;

    p = i * hl + j;
    index[kmer] = p;
  }
}

uint64_t apply_mask(const uint64_t &kmer, const std::string &mask) {
  uint64_t masked = 0;
  uint8_t base;
  uint8_t mp = 0;
  for (uint32_t i = 0; i < mask.size(); ++i) {
    if (mask[mask.size() - i - 1] == '1') {
      base = (kmer >> (2 * i)) & 3;
      masked |= (base << 2 * mp);
      ++mp;
    }
  }
  /**
  std::cout << kmer << " > " << masked << std::endl;
  std::cout << get_skmer(kmer, mask.size()) << " > " << get_skmer(masked, mp)
            << std::endl;
  **/
  return masked;
}

int main_kmc(int argc, char *argv[]) {
  // CLI
  int c;
  std::string mask = "";     // mask to use
  std::string out_dir = "."; // where to store all .npy
  bool with_mask = false;    // add mask to output file names
  opterr = 0;
  while ((c = getopt(argc, argv, "m:o:ph")) != -1) {
    switch (c) {
    case 'm':
      mask = optarg;
      break;
    case 'o':
      out_dir = optarg;
      break;
    case 'p':
      with_mask = true;
      break;
    case 'h':
      std::cerr << USAGE_MESSAGE << std::endl;
      exit(EXIT_SUCCESS);
    default:
      std::cerr << USAGE_MESSAGE << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  if (mask == "")
    with_mask = false;
  if (argc - optind < 1) {
    std::cerr << USAGE_MESSAGE << std::endl;
    exit(EXIT_FAILURE);
  }
  char *fpaths = argv[optind++];

  // Gwtting k from first file
  CKMCFile kmer_db;
  uint32_t mode, min_counter, pref_len, sign_len, min_c, counter, klen;
  uint64_t tot_kmers, max_c;
  std::string line;
  std::ifstream infile(fpaths);
  if (infile.is_open()) {
    getline(infile, line);
    if (!kmer_db.OpenForListing(line)) {
      std::cerr << "ERROR: cannot open " << line << std::endl;
      exit(1);
    }
    kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c,
                 tot_kmers);
    kmer_db.Close();
  } else {
    std::cerr << "ERROR: cannot open " << fpaths << std::endl;
    exit(1);
  }
  infile.close();

  if (klen > 15) {
    std::cerr << "ERROR: cannot work for k>15. Current k=" << klen << std::endl;
    exit(EXIT_FAILURE);
  }

  // Check mask
  if (mask.compare("") == 0) {
    mask = std::string(klen, '1');
  } else {
    if (mask.size() != klen) {
      std::cerr << "ERROR: mask and k error. Mask size: " << mask.size()
                << ", k: " << klen << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  uint32_t masked_klen = 0;
  for (const char &bit : mask)
    masked_klen += bit == '1';

  std::cerr << "Building index for k=" << klen << ". Mask: " << mask
            << ". Masked k=" << masked_klen << std::endl;

  // Init FCGR
  uint32_t fcgr_l = (1 << 2 * masked_klen); // total size of FCGR
  uint32_t hl = (1 << masked_klen);         // number of rows/columns
  std::vector<uint64_t> index(
      fcgr_l); // "map" from 2bit kmer to position in FCGR
  fill_index(index, masked_klen);

  uint64_t masked_kmer;
  std::vector<uint32_t> output(fcgr_l); // FCGR

  // Iterate over KMC databases
  CKmerAPI kmer_obj(klen);
  std::vector<uint64_t> kmer_l;
  infile.open(fpaths);
  std::filesystem::create_directories(out_dir);
  while (getline(infile, line)) {
    // Clear FCGR
    for (uint32_t i = 0; i < fcgr_l; ++i)
      output[i] = 0;

    std::cerr << "Parsing " << line << std::endl;
    if (!kmer_db.OpenForListing(line)) {
      std::cerr << "ERROR: cannot open " << line << std::endl;
      exit(1);
    }
    kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c,
                 tot_kmers);
    while (kmer_db.ReadNextKmer(kmer_obj, counter)) {
      kmer_obj.to_long(kmer_l);
      masked_kmer = apply_mask(kmer_l[0], mask);
      /**
      if (klen == masked_klen)
        assert(kmer_l[0] == masked_kmer);
      **/
      output[index[masked_kmer]] += counter;
    }
    std::string out_fn = out_dir + "/" +
                         std::filesystem::path(line).filename().string() +
                         (with_mask ? "." + mask : "") + ".npy";
    cnpy::npy_save(out_fn, &output[0], {hl, hl}, "w");
    std::fill(output.begin(), output.end(), 0);
    kmer_db.Close();
  }
  infile.close();

  return 0;
}

int main_fa(int argc, char *argv[]) {
  // CLI
  int c;
  uint8_t klen = 7;
  int n = -1;                // how many entries to consider from FASTA file
  std::string mask = "";     // mask to use
  std::string out_dir = "."; // where to store all .npy
  opterr = 0;
  while ((c = getopt(argc, argv, "k:n:m:o:h")) != -1) {
    switch (c) {
    case 'k':
      klen = std::stoi(optarg);
      break;
    case 'n':
      n = std::stoi(optarg);
      break;
    case 'm':
      mask = optarg;
      break;
    case 'o':
      out_dir = optarg;
      break;
    case 'h':
      // std::cerr << USAGE_MESSAGE << std::endl;
      exit(EXIT_SUCCESS);
    default:
      // std::cerr << USAGE_MESSAGE << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  if (argc - optind < 1) {
    std::cerr << USAGE_MESSAGE << std::endl;
    exit(EXIT_FAILURE);
  }
  char *fa_fn = argv[optind++];

  if (klen > 15) {
    std::cerr << "ERROR: cannot work for k>15. Current k=" << klen << std::endl;
    exit(EXIT_FAILURE);
  }

  // Check mask
  bool with_mask = true; // add mask to output file names
  if (mask.compare("") == 0) {
    with_mask = false;
    mask = std::string(klen, '1');
  } else {
    if (mask.size() != klen) {
      std::cerr << "ERROR: mask and k error. Mask size: " << mask.size()
                << ", k: " << klen << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  uint32_t masked_klen = 0;
  for (const char &bit : mask)
    masked_klen += bit == '1';

  std::cerr << "Building index for k=" << klen << ". Mask: " << mask
            << ". Masked k=" << masked_klen << std::endl;

  // Init FCGR
  uint32_t fcgr_l = (1 << 2 * masked_klen); // total size of FCGR
  uint32_t hl = (1 << masked_klen);         // number of rows/columns
  std::vector<uint64_t> index(
      fcgr_l); // "map" from 2bit kmer to position in FCGR
  fill_index(index, masked_klen);

  uint64_t masked_kmer;
  std::vector<uint32_t> output(fcgr_l); // FCGR
  for (uint32_t i = 0; i < fcgr_l; ++i)
    output[i] = 0;

  std::filesystem::create_directories(out_dir);

  // Iterate over FASTA
  gzFile fp = gzopen(fa_fn, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  int ee = 0;
  while ((l = kseq_read(seq)) >= 0) {
    if (n != -1 && ee == n)
      break;
    size_t pos = 0;      // current position
    char kmer[klen + 1]; // first kmer on sequence (plain)
    uint64_t kmer_d = 0; // kmer
    // uint64_t rckmer_d = 0; // reverse and complemented kmer
    uint64_t ckmer_d = 0; // canonical kmer
    uint8_t c;            // new character to append

    // first kmer
    strncpy(kmer, seq->seq.s, klen);
    kmer_d = k2d(kmer, klen);
    // rckmer_d = rc(kmer_d, klen);
    // ckmer_d = std::min(kmer_d, rckmer_d);
    ckmer_d = kmer_d;

    masked_kmer = apply_mask(ckmer_d, mask);
    output[index[masked_kmer]] += 1;

    // all other kmers
    ++pos;
    for (; pos < seq->seq.l - klen + 1; ++pos) {
      c = to_int[(uint8_t)seq->seq.s[pos]];
      kmer_d = lsappend(kmer_d, c, klen);
      // rckmer_d = rsprepend(rckmer_d, reverse_char(c), klen);
      // ckmer_d = std::min(kmer_d, rckmer_d);
      ckmer_d = kmer_d;

      masked_kmer = apply_mask(ckmer_d, mask);
      output[index[masked_kmer]] += 1;
    }

    // output
    std::string out_fn =
        out_dir + "/" + std::filesystem::path(seq->name.s).filename().string() +
        (with_mask ? "." + mask : "") + ".npy";
    cnpy::npy_save(out_fn, &output[0], {hl, hl}, "w");
    std::fill(output.begin(), output.end(), 0);

    ++ee;
  }
  kseq_destroy(seq);
  gzclose(fp);

  return 0;
}

int main(int argc, char *argv[]) {
  if (strcmp(argv[1], "fa") == 0)
    return main_fa(argc - 1, argv + 1);
  else
    return main_kmc(argc - 1, argv + 1);
  return 1;
}
