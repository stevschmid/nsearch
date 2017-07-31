#pragma once

#include "Sequence.h"
#include "TextReader.h"
#include "Utils.h"

#include <memory>

class SequenceReader {
public:
  SequenceReader( const std::string& pathToFile )
      : mTextReader( new TextFileReader( pathToFile ) ) {}

  SequenceReader( std::istream& is )
      : mTextReader( new TextStreamReader( is ) ) {}

  bool EndOfFile() const {
    return mTextReader->EndOfFile();
  }

  size_t NumBytesRead() const {
    return mTextReader->NumBytesRead();
  }

  size_t NumBytesTotal() const {
    return mTextReader->NumBytesTotal();
  }

  virtual void operator>>( Sequence& seq ) = 0;

  void Read( SequenceList& out, const size_t count ) {
    Sequence seq;

    for( size_t i = 0; i < count && !EndOfFile(); i++ ) {
      *this >> seq;
      out.push_back( std::move( seq ) );
    }
  }

protected:
  std::unique_ptr< TextReader > mTextReader;
};
