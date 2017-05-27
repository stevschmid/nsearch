#pragma once

#include "Sequence.h"
#include "Utils.h"

class LocalAligner {
public:
  LocalAligner( int matchScore = 1, int mismatchScore = -2,
      int gapOpenPenalty = 10, int gapExtendPenalty = 1 );
  void Align( const Sequence &query, const Sequence& target );

private:
  int mGapOpenPenalty, mGapExtendPenalty;
  int8_t mScoringMatrix[ NUC_MATRIX_SIZE * NUC_MATRIX_SIZE ];
};
