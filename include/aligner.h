#include <vector>
#include <string>
#include <cstdio>


#include <seqan/score.h>
#include <seqan/graph_msa.h>
#include <seqan/basic.h>


#include "aligner_opts.h"


#define MAXLINESZ 20480


// Seqan-compatible type
struct IntegerType_ {};
typedef seqan::SimpleType<uint16_t, IntegerType_> IntegerType;
typedef seqan::String<IntegerType> IntegerString;

struct IntegerSequence {
    std::string name;
    std::vector<uint16_t> data;
};

std::vector<IntegerSequence> read_integer_sequences_from_file(const std::string &input_fp);
seqan::StringSet<IntegerString> prepare_alignment_data(const std::vector<IntegerSequence> &integer_sequences);
seqan::Align<IntegerString> perform_alignment(seqan::StringSet<IntegerString> &integer_set);
void write_alignment(const seqan::Align<IntegerString> &alignment, const std::vector<IntegerSequence> &integer_sequences, const std::string &output_fp);
