#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "cnpy.h"
#include "kmc_api/kmc_file.h"

static const char *const USAGE_MESSAGE =
    "Usage: fcgr [-m MASK] <in-kmc-list>\n"
    "Options:\n"
    "        -m MASK       use this mask (default: 1{k})\n"
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

int main(int argc, char *argv[]) {
  // CLI
  int c;
  std::string mask = "";
  opterr = 0;
  while ((c = getopt(argc, argv, "m:h")) != -1) {
    switch (c) {
    case 'm':
      mask = optarg;
      break;
    case 'h':
      std::cerr << USAGE_MESSAGE << std::endl;
      exit(EXIT_SUCCESS);
    default:
      std::cerr << USAGE_MESSAGE << std::endl;
      exit(EXIT_FAILURE);
    }
  }
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

  uint32_t fcgr_l = (1 << 2 * masked_klen); // total size of FCGR
  uint32_t hl = (1 << masked_klen);         // number of rows/columns
  std::vector<uint64_t> index(
      fcgr_l); // "map" from 2bit kmer to position in FCGR
  fill_index(index, masked_klen);

  uint64_t masked_kmer;
  std::vector<uint32_t> output(fcgr_l); // FCGR

  CKmerAPI kmer_obj(klen);
  std::vector<uint64_t> kmer_l;
  infile.open(fpaths);
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
    cnpy::npy_save(line + ".npy", &output[0], {hl, hl}, "w");
    std::fill(output.begin(), output.end(), 0);
    kmer_db.Close();
  }
  infile.close();

  return 0;
}
