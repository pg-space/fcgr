#ifndef PS_KMER_H
#define PS_KMER_H

#include <stdint.h>
#include <stdlib.h>

static const uint8_t to_int[128] = {
    0, 0,       1, 2,       3,       0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,       0, 0,       0,       0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,       0, 0,       0,       0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,       0, 0,       0,       0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0,
    0, 0 /*A*/, 0, 1 /*C*/, 0,       0, 0, 2 /*G*/, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,       0, 0,       3 /*T*/, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0,
    0, 0 /*a*/, 0, 1 /*c*/, 0,       0, 0, 2 /*g*/, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,       0, 0,       3 /*t*/, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0};

uint8_t reverse_char(uint8_t c);

void d23(uint64_t kmer, int k, char *kk);

void d2s(uint64_t kmer, int k, char *kk);

uint64_t k2d(const char *kmer, uint8_t k);

uint64_t rc(uint64_t kmer, uint8_t k);

// left shift and append
uint64_t lsappend(uint64_t kmer, uint64_t c, uint8_t k);

// right shift and prepend
uint64_t rsprepend(uint64_t kmer, uint64_t c, uint8_t k);

#endif
