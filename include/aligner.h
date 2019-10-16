#include <vector>
#include <string>
#include <cstdio>


#include <seqan/score.h>
#include <seqan/graph_msa.h>
#include <seqan/basic.h>


#include "aligner_opts.h"


#define MAXLINESZ 20480


// Seqan-compatible types
struct IntegerType_ {};
typedef seqan::SimpleType<uint16_t, IntegerType_> IntegerType;
typedef seqan::String<IntegerType> IntegerString;
typedef seqan::StringSet<IntegerString> IntegerStringSet;
typedef seqan::StringSet<IntegerString, seqan::Dependent<>> DepIntegerStringSet;
typedef seqan::Graph<seqan::Alignment<DepIntegerStringSet>> IntegerAlignGraph;

struct IntegerSequence {
    std::string name;
    std::vector<uint16_t> data;
};

struct GraphData {
    IntegerAlignGraph graph;
    IntegerStringSet sequence_set;
};

std::vector<IntegerSequence> read_integer_sequences_from_file(const std::string &input_fp);
seqan::StringSet<IntegerString> prepare_alignment_data(const std::vector<IntegerSequence> &integer_sequences);
GraphData perform_alignment(IntegerStringSet &integer_set, const arguments::ScoringScheme score);
void write_alignment(const IntegerAlignGraph &graph, const std::vector<IntegerSequence> integer_sequences, const std::string &output_fp);
