#include "aligner.h"


int main(int argc, char *argv[]) {
    // Get command line arguments
    arguments::Arguments args = arguments::get_arguments(argc, argv);
    // Read integer sequences from file and prepare alignment data
    IntegerSequences integer_sequences = read_integer_sequences_from_file(args.input_fp);
    seqan::StringSet<IntegerString> integer_set = prepare_alignment_data(integer_sequences.sequences);
    // Perform alignment and write out
    GraphData graph_data = perform_alignment(integer_set, args.scoring_scheme, args.use_global);
    write_alignment(graph_data.graph, integer_sequences, args.output_fp);
    return 0;
}

IntegerSequences read_integer_sequences_from_file(const std::string &input_fp) {
    // TODO: better file handling, catch literally any exception instead of allowing prog to segfault
    // Read data into structs
    FILE *input_fh = fopen(input_fp.c_str(), "r");
    char line[MAXLINESZ];
    size_t position = 0;
    IntegerSequences integer_sequences;
    while(fgets(line, MAXLINESZ, input_fh)) {
        // First pop sequence name from line
        IntegerSequence integer_sequence;
        std::string line_string(line), token;
        position = line_string.find(",");
        integer_sequence.name = line_string.substr(0, position);
        line_string.erase(0, position + 1);
        integer_sequences.sequences.push_back(integer_sequence);
        // Get each token and process
        while ((position = line_string.find(",")) != std::string::npos) {
            token = line_string.substr(0, position);
            line_string.erase(0, position + 1);
            process_token(token, integer_sequences);
        }
        // Process last token
        process_token(line_string, integer_sequences);
    }
    fclose(input_fh);
    return integer_sequences;
}

void process_token(const std::string token, IntegerSequences &integer_sequences) {
    uint16_t id_num;
    uint16_t value = std::stoi(token);
    if (integer_sequences.value_map.find(value) == integer_sequences.value_map.end()) {
        if (integer_sequences.id_map.size() == IDMAX) {
            fprintf(stderr, "error: exceeded hard coded limit of %d integers\n", IDMAX);
            exit(1);
        }
        id_num = integer_sequences.id_map.size();
        if (id_num == 45) {
            id_num = 0xFF;
        }
        integer_sequences.value_map.emplace(value, id_num);
        integer_sequences.id_map.emplace(id_num, value);
    } else {
        id_num = integer_sequences.value_map.find(value)->second;
    }
    integer_sequences.sequences.back().data.push_back(id_num);
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

GraphData perform_alignment(IntegerStringSet &integer_set, const arguments::ScoringScheme &score, bool use_global) {
    // Place graph data into struct, stringSet must stay in scope
    GraphData graph_data = {IntegerAlignGraph(integer_set), seqan::stringSet(graph_data.graph)};
    // Create name array
    seqan::String<seqan::String<char>> names;
    seqan::resize(names, seqan::length(graph_data.sequence_set), seqan::String<char>("empty"));
    // Set alignment options
    seqan::MsaOptions<IntegerString, seqan::Score<int>> msa_options;
    msa_options.sc = seqan::Score<int>(score.match, score.mismatch, score.gapextend, score.gapopen);
    if (use_global) {
        // Set global alignment
        seqan::appendValue(msa_options.method, 0);
    }
    // Set longest subsequence
    seqan::appendValue(msa_options.method, 3);
    // Run alignment and return
    seqan::globalMsaAlignment(graph_data.graph, graph_data.sequence_set, names, msa_options);
    return graph_data;
}

void write_alignment(const IntegerAlignGraph &graph, const IntegerSequences integer_sequences, const std::string &output_fp) {
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
        fprintf(output_fh, integer_sequences.sequences[row].name.c_str());
        for (TSize col = 0; col < cols; ++col) {
            uint16_t id_num = alignment[row * cols + col];
            // Print gaps or values; ascii 45 == '-'
            if (id_num == 45) {
                fprintf(output_fh, "\t-");
            } else {
                uint16_t value = integer_sequences.id_map.at(id_num);
                if (value == 0xFF) {
                    value = 45;
                }
                fprintf(output_fh, "\t%d", value);
            }
        }
        fprintf(output_fh, "\n");
    }
    fclose(output_fh);
}
