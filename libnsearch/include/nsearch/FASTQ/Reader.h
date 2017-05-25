#pragma once

#include "../TextReader.h"
#include "../Sequence.h"

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
      std::string identifier, sequence, plusline, quality;

      (*mTextReader) >> identifier;
      (*mTextReader) >> sequence;
      (*mTextReader) >> plusline;
      (*mTextReader) >> quality;

      seq = Sequence( identifier.substr( 1 ), sequence, quality );
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
