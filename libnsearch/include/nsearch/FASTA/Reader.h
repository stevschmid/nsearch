#pragma once

#include "../SequenceReader.h"

namespace FASTA {
  class Reader : public SequenceReader {
  public:
    using SequenceReader::SequenceReader;

    void operator>>( Sequence &seq ) {
      std::string identifier, sequence;
      if( mLastLine.empty() ) {
        (*mTextReader) >> identifier;
      } else {
        identifier = mLastLine;
      }

      std::string line;
      while( !EndOfFile() ) {
        (*mTextReader) >> line;

        if( line[ 0 ] == '>' ) {
          mLastLine = line;
          break;
        }

        sequence += line;
      }

      UpcaseString( sequence );
      seq = Sequence( identifier.substr( 1 ), sequence );
    }

  private:
    std::string mLastLine;
  };
}
