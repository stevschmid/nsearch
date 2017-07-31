#pragma once

#include "../Sequence.h"

#include <fstream>

namespace FASTQ {
class Writer {
private:
  std::ofstream mFile;
  std::ostream& mOutput;

public:
  Writer( std::ostream& output ) : mOutput( output ) {}

  Writer( const std::string& pathToFile )
      : mFile( pathToFile ), mOutput( mFile ) {}

  void operator<<( const Sequence& seq ) {
    mOutput << '@' << seq.identifier << std::endl;
    mOutput << seq.sequence << std::endl;
    mOutput << '+' << std::endl;
    mOutput << seq.quality << std::endl;
  }
};
} // namespace FASTQ
