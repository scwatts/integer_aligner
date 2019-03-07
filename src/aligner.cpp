#include <vector>
#include <string>
#include <cstdio>


#include <seqan/graph_msa.h>
#include <seqan/basic.h>


#include "aligner_opts.h"


#define MAXLINESZ 20480


// Define seqan-compatible type
struct IntegerType_ {};
typedef seqan::SimpleType<uint16_t, IntegerType_> IntegerType;
typedef seqan::String<IntegerType> IntegerString;

struct IntegerSequence {
    std::string name;
    std::vector<uint16_t> data;
};

// Function declarations
std::vector<IntegerSequence> read_integer_sequences_from_file(const std::string &input_fp);
seqan::StringSet<IntegerString> prepare_alignment_data(const std::vector<IntegerSequence> &integer_sequences);
seqan::Align<IntegerString> perform_alignment(seqan::StringSet<IntegerString> &integer_set);
void write_alignment(const seqan::Align<IntegerString> &alignment, const std::vector<IntegerSequence> &integer_sequences, const std::string &output_fp);


int main(int argc, char *argv[]) {
    // Get command line arguments
    arguments::Arguments args = arguments::get_arguments(argc, argv);
    // Read integer sequences from file and prepare alignment data
    std::vector<IntegerSequence> integer_sequences = read_integer_sequences_from_file(args.input_fp);
    seqan::StringSet<IntegerString> integer_set = prepare_alignment_data(integer_sequences);
    // Perform alignment and write out
    seqan::Align<IntegerString> alignment = perform_alignment(integer_set);
    write_alignment(alignment, integer_sequences, args.output_fp);
    return 0;
}

std::vector<IntegerSequence> read_integer_sequences_from_file(const std::string &input_fp) {
    // TODO: better file handling, catch literally any exception instead of allowing prog to segfault
    // Declare return value
    std::vector<IntegerSequence> integer_sequences;
    // Read data into structs
    FILE *input_fh = fopen(input_fp.c_str(), "r");
    char line[MAXLINESZ];
    size_t position = 0;
    while(fgets(line, MAXLINESZ, input_fh)) {
        // First pop sequence name from line
        IntegerSequence integer_sequence;
        std::string line_string(line), token;
        position = line_string.find(",");
        integer_sequence.name = line_string.substr(0, position);
        line_string.erase(0, position + 1);
        while ((position = line_string.find(",")) != std::string::npos) {
            // Convert string token to integer and clip from line string
            token = line_string.substr(0, position);
            integer_sequence.data.push_back(std::stoi(token));
            line_string.erase(0, position + 1);
        }
        // Add last element
        integer_sequence.data.push_back(std::stoi(line_string));
        // Update return vector with new integer sequence
        integer_sequences.push_back(integer_sequence);
    }
    fclose(input_fh);
    return integer_sequences;
}

seqan::StringSet<IntegerString> prepare_alignment_data(const std::vector<IntegerSequence> &integer_sequences) {
    // TODO: there is many copies here - could be more efficient
    seqan::StringSet<IntegerString> integer_set;
    for (const auto& integer_sequence : integer_sequences) {
        IntegerString integer_string;
        for (const auto& v : integer_sequence.data) {
            seqan::appendValue(integer_string, v);
        }
        seqan::appendValue(integer_set, integer_string);
    }
    return integer_set;
}

seqan::Align<IntegerString> perform_alignment(seqan::StringSet<IntegerString> &integer_set) {
    seqan::Align<IntegerString> alignment(integer_set);
    seqan::globalMsaAlignment(alignment, seqan::EditDistanceScore());
    return alignment;
}

void write_alignment(const seqan::Align<IntegerString> &alignment, const std::vector<IntegerSequence> &integer_sequences, const std::string &output_fp) {
    // Needed for readable code
    typedef seqan::Align<IntegerString> const TAlign;
    typedef typename seqan::Position<typename seqan::Rows<TAlign>::Type>::Type TRowsPosition;
    typedef typename seqan::Position<TAlign>::Type TPosition;
    typedef typename seqan::Iterator<typename seqan::Row<TAlign>::Type const, seqan::Standard>::Type TIter;
    // Iterate and write result
    FILE *output_fh = fopen(output_fp.c_str(), "w");
    TRowsPosition row_count = seqan::length(seqan::rows(alignment));
    TPosition end_position = std::min(seqan::length(seqan::row(alignment, 0)), seqan::length(seqan::row(alignment, 1)));
    // Iterate rows
    for (long unsigned int i = 0; i < row_count; ++i) {
        // Print sequence name
        fprintf(output_fh, "%s", integer_sequences[i].name.c_str());
        // Iterate integers of row
        auto& row_data = seqan::row(alignment, i);
        TIter row_iter = iter(row_data, 0);
        for (int j = 0; row_iter != iter(row_data, end_position); ++row_iter, ++j) {
            if (seqan::isGap(row_iter)) {
                fprintf(output_fh, "\t-");
            } else {
                // TODO: there's probably a better way than this ¯\_(ツ)_/¯
                uint16_t value = (row_iter._container->_source.data_value->data_begin+row_iter._sourcePosition)->value;
                fprintf(output_fh, "\t%d", value);
            }
        }
        fprintf(output_fh, "\n");
    }
    fclose(output_fh);
}
