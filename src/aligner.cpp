#include "aligner.h"


#define OFFSET 0xFF


int main(int argc, char *argv[]) {
    // Get command line arguments
    arguments::Arguments args = arguments::get_arguments(argc, argv);
    // Read integer sequences from file and prepare alignment data
    std::vector<IntegerSequence> integer_sequences = read_integer_sequences_from_file(args.input_fp);
    seqan::StringSet<IntegerString> integer_set = prepare_alignment_data(integer_sequences);
    // Perform alignment and write out
    GraphData graph_data = perform_alignment(integer_set, args.scoring_scheme);
    write_alignment(graph_data.graph, integer_sequences, args.output_fp);
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
            line_string.erase(0, position + 1);
            // Convert and shift into non-ascii range
            integer_sequence.data.push_back(std::stoi(token) + OFFSET);
        }
        integer_sequence.data.push_back(std::stoi(line_string) + OFFSET);
        // Update return vector with new integer sequence
        integer_sequences.push_back(integer_sequence);
    }
    fclose(input_fh);
    return integer_sequences;
}

seqan::StringSet<IntegerString> prepare_alignment_data(const std::vector<IntegerSequence> &integer_sequences) {
    // TODO: there are many copies here - could be more efficient
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

GraphData perform_alignment(IntegerStringSet &integer_set, const arguments::ScoringScheme score) {
    // Place graph data into struct, stringSet must stay in scope
    GraphData graph_data = {IntegerAlignGraph(integer_set), seqan::stringSet(graph_data.graph)};
    // Create name array
    seqan::String<seqan::String<char>> names;
    seqan::resize(names, seqan::length(graph_data.sequence_set), seqan::String<char>("empty"));
    // Set options - run pairwise global and local, and overlap, and longest subsequence
    seqan::MsaOptions<IntegerString, seqan::Score<int>> msa_options;
    msa_options.sc = seqan::Score<int>(score.match, score.mismatch, score.gapextend, score.gapopen);
    seqan::appendValue(msa_options.method, 0);  // Global pairwise
    // TODO: explicitly set alphabet size to use local pairwise
    //seqan::appendValue(msa_options.method, 1);  // Local pairwise
    //seqan::appendValue(msa_options.method, 2);  // Overlap
    //seqan::appendValue(msa_options.method, 3);  // Longest subsequence
    // Align
    seqan::globalMsaAlignment(graph_data.graph, graph_data.sequence_set, names, msa_options);
    return graph_data;
}

void write_alignment(const IntegerAlignGraph &graph, const std::vector<IntegerSequence> integer_sequences, const std::string &output_fp) {
    // Typedefs
    typedef typename seqan::Size<IntegerAlignGraph>::Type TSize;
    // Create alignment matrix from graphical alignment
    IntegerString alignment;
    seqan::convertAlignment(graph, alignment);
    // Write sequences
    FILE *output_fh = fopen(output_fp.c_str(), "w");
    TSize rows = seqan::length(seqan::stringSet(graph));
    TSize cols = seqan::length(alignment) / rows;
    for (TSize row = 0; row < rows; ++row) {
        fprintf(output_fh, integer_sequences[row].name.c_str());
        for (TSize col = 0; col < cols; ++col) {
            uint16_t value = alignment[row * cols + col];
            // Print gaps or values; ascii 45 == '-'
            if (value == 45) {
                fprintf(output_fh, "\t-");
            } else {
                fprintf(output_fh, "\t%d", value - OFFSET);
            }
        }
        fprintf(output_fh, "\n");
    }
    fclose(output_fh);
}
