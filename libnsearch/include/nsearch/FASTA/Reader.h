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
      static std::string lastLine;

      std::string identifier, sequence;
      if( lastLine.empty() ) {
        (*mTextReader) >> identifier;
      } else {
        identifier = lastLine;
      }

      std::string line;
      while( !EndOfFile() ) {
        (*mTextReader) >> line;

        if( line[ 0 ] == '>' ) {
          lastLine = line;
          break;
        }

        sequence += line;
      }

      seq = Sequence( identifier.substr( 1 ), toupper( sequence ) );
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
