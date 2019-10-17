#ifndef __ALIGNER_OPTS_H__
#define __ALIGNER_OPTS_H__


#include <cstdio>
#include <vector>
#include <string>


#include <getopt.h>


#include "config.h"  // autoconf header


namespace arguments {


struct ScoringScheme {
    int match = 0;
    int mismatch = -2;
    int gapextend = -1;
    int gapopen = -2;
};

struct Arguments {
    std::string input_fp;
    std::string output_fp;
    ScoringScheme scoring_scheme;
    bool use_global = false;
};

void print_help(FILE *stdst);
void print_version(FILE *stdst);
Arguments get_arguments(int argc, char **argv);
int int_from_optarg(const char *optarg);


} // namespace arguments


#endif
