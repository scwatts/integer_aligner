#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>


#include <seqan/score.h>
#include <seqan/graph_msa.h>
#include <seqan/basic.h>


#include "aligner_opts.h"


#define MAXLINESZ 20480
#define IDOFFSET 0xFF
#define IDMAX 0xFF


// Seqan-compatible types
struct IntegerType_ {};
typedef seqan::SimpleType<uint16_t, IntegerType_> IntegerType;
typedef seqan::String<IntegerType> IntegerString;
typedef seqan::StringSet<IntegerString> IntegerStringSet;
typedef seqan::StringSet<IntegerString, seqan::Dependent<>> DepIntegerStringSet;
typedef seqan::Graph<seqan::Alignment<DepIntegerStringSet>> IntegerAlignGraph;
// Other types
typedef std::unordered_map<uint16_t,uint16_t> uint16_map;

// Set alphabet size
template <>
struct seqan::ValueSize<IntegerType>
{
    typedef uint8_t Type;
    static const Type VALUE = IDMAX;
};

template <>
struct seqan::BitsPerValue<IntegerType>
{
    typedef uint8_t Type;
    static const uint8_t VALUE = sizeof(IDMAX);
};


// Data structs
struct IntegerSequence {
    std::string name;
    std::vector<uint16_t> data;
};

struct IntegerSequences {
    std::vector<IntegerSequence> sequences;
    uint16_map value_map;
    uint16_map id_map;
};

struct GraphData {
    IntegerAlignGraph graph;
    IntegerStringSet sequence_set;
};

// Functions
IntegerSequences read_integer_sequences_from_file(const std::string &input_fp);
void process_token(const std::string token, IntegerSequences &integer_sequences);
seqan::StringSet<IntegerString> prepare_alignment_data(const std::vector<IntegerSequence> &integer_sequences);
GraphData perform_alignment(IntegerStringSet &integer_set, const arguments::ScoringScheme &score, bool use_global);
void write_alignment(const IntegerAlignGraph &graph, const IntegerSequences integer_sequences, const std::string &output_fp);
