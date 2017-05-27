#pragma once

#include "../TextReader.h"
#include "../Sequence.h"
#include "../Utils.h"

#include <memory>

namespace FASTA {
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

    size_t NumBytesRead() const {
      return mTextReader->NumBytesRead();
    }

    size_t NumBytesTotal() const {
      return mTextReader->NumBytesTotal();
    }

  private:
    std::unique_ptr< TextReader > mTextReader;
    std::string mLastLine;
  };
}
