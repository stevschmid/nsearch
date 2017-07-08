#pragma once

#include "../Utils.h"

#define MAXINT INT_MAX/2 //prevent overflow
#define MININT -INT_MIN/2 //prevent underflow

enum class AlignmentDirection {
  forwards, backwards
};

enum class CigarOp : char {
  UNKNOWN   = ' ',
  MATCH     = 'M',
  MISMATCH  = 'X',
  DELETION  = 'D',
  INSERTION = 'I',
};

class CigarEntry {
public:
  int count = 0;
  CigarOp op = CigarOp::UNKNOWN;

  CigarEntry() { }
  CigarEntry( int count, CigarOp op )
    : count( count ), op( op ) { }
};

using Cigar = std::deque< CigarEntry >;
