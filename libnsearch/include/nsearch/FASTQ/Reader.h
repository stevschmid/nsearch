#pragma once

#include "../SequenceReader.h"

namespace FASTQ {
class Reader : public SequenceReader {
public:
  using SequenceReader::SequenceReader;

  void operator>>( Sequence& seq ) {
    ( *mTextReader ) >> seq.identifier;
    ( *mTextReader ) >> seq.sequence;
    ( *mTextReader ) >> seq.quality; // skip plusline
    ( *mTextReader ) >> seq.quality;

    // delete '>'
    seq.identifier.erase( seq.identifier.begin(), seq.identifier.begin() + 1 );

    UpcaseString( seq.sequence ); // atc -> ATC
    UpcaseString( seq.quality );
  }
};
} // namespace FASTQ
