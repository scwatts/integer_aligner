#include "aligner_opts.h"


namespace arguments {


void print_help(FILE *stdst) {
    // Print version
    print_version(stdst);
    // Print Usage
    fprintf(stdst, "\n");
    fprintf(stdst, "Usage:\n");
    fprintf(stdst, "  integer_aligner [options]\n");
    fprintf(stdst, "\n");
    fprintf(stdst, "IO:\n");
    fprintf(stdst, "    -i <filepath>, --input_fp <filepath>\n");
    fprintf(stdst, "                    Input alignment file containing integers to align\n");
    fprintf(stdst, "    -o <filepath>, --output_fp <filepath>\n");
    fprintf(stdst, "                    output filepath for aligned integers\n");
    fprintf(stdst, "\n");
    fprintf(stdst, "Alignment:\n");
    fprintf(stdst, "    -g, --global\n");
    fprintf(stdst, "                    Use global alignment (default: false)\n");
    fprintf(stdst, "\n");
    fprintf(stdst, "Scoring:\n");
    fprintf(stdst, "    -a <int>, --match <int>\n");
    fprintf(stdst, "                    Match score (default: 0)\n");
    fprintf(stdst, "    -b <int>, --mismatch <int>\n");
    fprintf(stdst, "                    Mismatch score (default: -2)\n");
    fprintf(stdst, "    -c <int>, --gapextend <int>\n");
    fprintf(stdst, "                    Gap extend score (default: -1)\n");
    fprintf(stdst, "    -d <int>, --gapopen <int>\n");
    fprintf(stdst, "                    Gap open score (default: -2)\n");
    fprintf(stdst, "\n");
    fprintf(stdst, "Other:\n");
    fprintf(stdst, "  -h, --help\n");
    fprintf(stdst, "                    display this help and exit\n");
    fprintf(stdst, "  -v, --version\n");
    fprintf(stdst, "                    display version information and exit\n");
}

void print_version(FILE *stdst) {
    fprintf(stdst, "Program: integer_aligner\n");
    fprintf(stdst, "Version %s\n", PACKAGE_VERSION);
    fprintf(stdst, "Contact: Stephen Watts (s.watts2@student.unimelb.edu.au)\n");
}

Arguments get_arguments(int argc, char **argv) {
    Arguments args;
    // Long opts
    struct option long_options[] =
        {
            {"input_fp", required_argument, NULL, 'i'},
            {"output_fp", required_argument, NULL, 'o'},
            {"global", no_argument, NULL, 'g'},
            {"match", optional_argument, NULL, 'a'},
            {"mismatch", optional_argument, NULL, 'b'},
            {"gapextend", optional_argument, NULL, 'c'},
            {"gapopen", optional_argument, NULL, 'd'},
            {"output_fp", required_argument, NULL, 'o'},
            {"version", no_argument, NULL, 'v'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, 0, 0}
        };

    // Parse commandline arguments
    while (1) {
        // Parser variables
        int long_options_index = 0;
        int c;

        // Parser
        c = getopt_long(argc, argv, "i:o:ga:b:c:d:vh", long_options, &long_options_index);

        // If no more arguments to parse, break
        if (c == -1) {
            break;
        }

        // Process current arguments
        switch(c) {
            case 'i':
                args.input_fp = optarg;
                break;
            case 'o':
                args.output_fp = optarg;
                break;
            case 'g':
                args.use_global = true;
                break;
            case 'a':
                args.scoring_scheme.match = int_from_optarg(optarg);
                break;
            case 'b':
                args.scoring_scheme.mismatch = int_from_optarg(optarg);
                break;
            case 'c':
                args.scoring_scheme.gapextend = int_from_optarg(optarg);
                break;
            case 'd':
                args.scoring_scheme.gapopen = int_from_optarg(optarg);
                break;
            case 'v':
                print_version(stdout);
                exit(0);
            case 'h':
                print_help(stdout);
                exit(0);
            default:
                exit(1);
        }
    }

    // Ensure that input filepaths are set
    if (args.input_fp.empty()) {
        print_help(stderr);
        fprintf(stderr,"\n%s: error: argument -i/--input_fp is required\n", argv[0]);
        exit(1);
    }
    if (args.output_fp.empty()) {
        print_help(stderr);
        fprintf(stderr,"\n%s: error: argument -o/--output_fp is required\n", argv[0]);
        exit(1);
    }

    return args;
}


int int_from_optarg(const char *optarg) {
    // Check at most the first 8 characters are numerical
    std::string optstring(optarg);
    int offset = (optstring.at(0) == '-') ? 1 : 0;
    std::string string_int = optstring.substr(0, 8-offset);
    for (std::string::iterator it = string_int.begin()+offset; it != string_int.end(); ++it) {
        if (!isdigit(*it)) {
            fprintf(stderr, "This doesn't look like a usable integer: %s\n", optarg);
            exit(1);
        }
    }
    return std::atoi(string_int.c_str());
}


} // namespace arguments
