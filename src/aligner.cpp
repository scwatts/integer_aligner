#include <seqan/graph_msa.h>
#include <seqan/basic.h>


#include "aligner_opts.h"



int main(int argc, char *argv[]) {
    // Get command line arguments
    arguments::Arguments args = arguments::get_arguments(argc, argv);


    // Prepare integer sequences for alignment
    struct Int_ {};
    typedef seqan::SimpleType<uint16_t, Int_> Int;
    typedef seqan::String<Int> IntString;

    std::vector<uint16_t> val_a = {1,8,7,7,7,7,7,7,1000,4,5};
    std::vector<uint16_t> val_b = {2,8,9,9,4,5};
    std::vector<uint16_t> val_c = {3,8,9,9,4,5};
    IntString test_a, test_b, test_c;
    for (const auto v : val_a) {
        seqan::appendValue(test_a, v);
    }
    for (const auto v : val_b) {
        seqan::appendValue(test_b, v);
    }
    for (const auto v : val_c) {
        seqan::appendValue(test_c, v);
    }

    seqan::StringSet<IntString> intSet;
    seqan::appendValue(intSet, test_a);
    seqan::appendValue(intSet, test_b);
    seqan::appendValue(intSet, test_c);

    // Perform alignment
    seqan::Align<IntString> align_int(intSet);
    seqan::globalMsaAlignment(align_int, seqan::EditDistanceScore());

    // Write out alignment
    typedef seqan::Align<IntString> const TAlign;
    typedef typename seqan::Position<typename seqan::Rows<TAlign>::Type>::Type TRowsPosition;
    typedef typename seqan::Position<TAlign>::Type TPosition;

    TRowsPosition row_count = seqan::length(seqan::rows(align_int));
    TPosition end_ = std::min(seqan::length(seqan::row(align_int, 0)), seqan::length(seqan::row(align_int, 1)));

    for (long unsigned int i = 0; i < 2 * row_count - 1; ++i) {
        // Print integers
        if ((i % 2) == 0) {
            auto & row_ = seqan::row(align_int, i / 2);
            typedef typename seqan::Iterator<typename seqan::Row<TAlign>::Type const, seqan::Standard>::Type TIter;
            TIter begin1_ = iter(row_, 0);
            TIter end1_ = iter(row_, end_);
            for (; begin1_ != end1_; ++begin1_) {
                if (begin1_._sourcePosition > 0) fprintf(stdout, "\t");
                if (seqan::isGap(begin1_)) {
                    fprintf(stdout, "-");
                } else {
                    uint16_t value = (begin1_._container->_source.data_value->data_begin+begin1_._sourcePosition)->value;
                    fprintf(stdout, "%d", value);
                }
            }
        // Aligned status (match, gap, non-match)
        } else {
            for (unsigned int j = 0; j < end_; ++j) {
                if (j > 0) fprintf(stdout, "\t");
                if ((!seqan::isGap(seqan::row(align_int, (i - 1) / 2), j)) &&
                    (!seqan::isGap(seqan::row(align_int, (i + 1) / 2), j)) &&
                    (seqan::row(align_int, (i - 1) / 2)[j] == seqan::row(align_int, (i + 1) / 2)[j])) {
                    fprintf(stdout, "|");
                } else {
                    fprintf(stdout, " ");
                }
            }
        }
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");

    return 0;
}
