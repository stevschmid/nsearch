#pragma once

#include "../TextReader.h"
#include "../Sequence.h"
#include "../Utils.h"

#include <memory>

namespace FASTQ {
  class Reader {
  public:
    Reader( const std::string &pathToFile )
      : mTextReader( new TextFileReader( pathToFile ) )
    {
    }

    Reader( std::istream &is )
      : mTextReader( new TextStreamReader( is ) )
    {
    }

    bool EndOfFile() const {
      return mTextReader->EndOfFile();
    }

    void operator>>( Sequence &seq ) {
      (*mTextReader) >> seq.identifier;
      (*mTextReader) >> seq.sequence;
      (*mTextReader) >> seq.quality; // skip plusline
      (*mTextReader) >> seq.quality;

      // delete '>'
      seq.identifier.erase( seq.identifier.begin(), seq.identifier.begin() + 1 );

      UpcaseString( seq.sequence ); // atc -> ATC
      UpcaseString( seq.quality );
    }

    size_t NumBytesRead() const {
      return mTextReader->NumBytesRead();
    }

    size_t NumBytesTotal() const {
      return mTextReader->NumBytesTotal();
    }

  private:
    std::unique_ptr< TextReader > mTextReader;
  };
}
