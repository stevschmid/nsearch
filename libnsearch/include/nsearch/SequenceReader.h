#pragma once

#include "TextReader.h"
#include "Sequence.h"
#include "Utils.h"

#include <memory>

class SequenceReader {
public:
  SequenceReader( const std::string &pathToFile )
    : mTextReader( new TextFileReader( pathToFile ) )
  {
  }

  SequenceReader( std::istream &is )
    : mTextReader( new TextStreamReader( is ) )
  {
  }

  bool EndOfFile() const {
    return mTextReader->EndOfFile();
  }

  size_t NumBytesRead() const {
    return mTextReader->NumBytesRead();
  }

  size_t NumBytesTotal() const {
    return mTextReader->NumBytesTotal();
  }

  virtual void operator>>( Sequence &seq ) = 0;

protected:
  std::unique_ptr< TextReader > mTextReader;
};
