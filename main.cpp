#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>

#include "cnpy.h"
#include "kmc_api/kmc_file.h"

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

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "ERROR. Usage: fcgr <in-kmc-list>" << std::endl;
    exit(1);
  }
 
  char *fpaths = argv[1];

  CKMCFile kmer_db;
  uint32_t mode, min_counter, pref_len, sign_len, min_c, counter, klen;
  uint64_t tot_kmers, max_c;

  std::string line;
  std::ifstream infile(fpaths);
  if (infile.is_open()) {
    std::cerr << "Get k from first file" << std::endl;
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

  std::cerr << "Building index for k=" << klen << std::endl;

  uint32_t l = (1 << 2 * klen);
  uint32_t hl = (1 << klen);
  std::vector<uint64_t> index(l);
  fill_index(index, klen);

  std::vector<uint32_t> output(l);
  CKmerAPI kmer_obj(klen);
  std::vector<uint64_t> kmer_l;
  infile.open(fpaths);
  while (getline(infile, line)) {
    std::cerr << "Parsing " << line << std::endl;
    if (!kmer_db.OpenForListing(line)) {
      std::cerr << "ERROR: cannot open " << line << std::endl;
      exit(1);
    }
    kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c,
                 tot_kmers);
    while (kmer_db.ReadNextKmer(kmer_obj, counter)) {
      kmer_obj.to_long(kmer_l);
      output[index[kmer_l[0]]] = counter;
    }
    cnpy::npy_save(line + ".npy", &output[0], {hl, hl}, "w");
    std::fill(output.begin(), output.end(), 0);
    kmer_db.Close();
  }
  infile.close();

  return 0;
}
