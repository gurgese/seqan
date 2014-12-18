#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
    BamFileIn bamFileIn("example.bam");

    BamHeader header;
    readRecord(header, bamFileIn);

    typedef SmartFileContext<BamFileIn, void>::Type TBamContext;

    TBamContext const & bamContext = context(bamFileIn);

    for (unsigned i = 0; i < length(contigNames(bamContext)); ++i)
        std::cout << contigNames(bamContext)[i] << '\t'
                  << contigLengths(bamContext)[i] << '\n';

    return 0;
}
