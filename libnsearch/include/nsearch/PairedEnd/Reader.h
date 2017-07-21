#pragma once

#include "../Sequence.h"
#include "../FASTQ/Reader.h"

namespace PairedEnd {
  class Reader {
  private:
    FASTQ::Reader mFwdReader;
    FASTQ::Reader mRevReader;

  public:
    Reader( const std::string &pathToFwdFile, const std::string &pathToRevFile )
      : mFwdReader( pathToFwdFile ), mRevReader( pathToRevFile )
    {
    }

    Reader( std::istream &fwd, std::istream &rev )
      : mFwdReader( fwd ), mRevReader( rev )
    {
    }

    void Read( SequenceList &fwds, SequenceList &revs, int count ) {
      for( int i = 0; i < count; i++ ) {
        Sequence fwd, rev;

        if( !Read( fwd, rev ) )
          break;

        fwds.push_back( fwd );
        revs.push_back( rev );
      }
    }

    bool Read( Sequence &fwd, Sequence &rev ) {
      if( EndOfFile() )
        return false;

      mFwdReader >> fwd;
      mRevReader >> rev;

      return true;
    }

    bool EndOfFile() const {
      return mFwdReader.EndOfFile() || mRevReader.EndOfFile();
    }

    std::streampos NumBytesRead() const {
      return mFwdReader.NumBytesRead() + mRevReader.NumBytesRead();
    }

    std::streampos NumBytesTotal() const {
      return mFwdReader.NumBytesTotal() + mRevReader.NumBytesTotal();
    }
  };
}
