#ifndef PALSS_USAGE_H
#define PALSS_USAGE_H

static const char *const VERSION = "fcgr, v1.0.0";

static const char *const USAGE_MESSAGE = "Usage: fcgr [single|multi] -h\n";

static const char *const SINGLE_USAGE_MESSAGE =
    "Usage: fcgr single <in-kmc-list>\n"
    "Options:\n"
    "        -m MASK       use this mask (default: 1{k})\n"
    "        -o DIR        store output(s) in this directory (default: ./)\n"
    "        -p            add mask to output filenames (default: false)\n"
    "        -h            display this help and exit\n";

static const char *const MULTI_USAGE_MESSAGE =
    "Usage: fcgr multi <in.fasta>\n"
    "Options:\n"
    "        -m MASK       use this mask (default: 1{k})\n"
    "        -n INT        how many entries to consider from input FASTA file "
    "(default: -1, all entries)\n"
    "        -o DIR        store output(s) in this directory (default: ./)\n"
    "        -p            add mask to output filenames (default: false)\n"
    "        -h            display this help and exit\n";
#endif
