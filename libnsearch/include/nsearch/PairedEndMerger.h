#pragma once

#include "Sequence.h"

class PairedEndMerger
{
public:
  PairedEndMerger( int minOverlap, double minIdentity );
  bool Merge( Sequence& merged, const Sequence &fwd, const Sequence &rev ) const;

private:
  int mMinOverlap;
  double mMinIdentity;

  typedef struct {
    size_t length;
    size_t pos1;
    size_t pos2;
  } overlapInfo;

  int ComputeOverlapScore( const char *sequence1, const char *sequence2, size_t len ) const;
  bool FindBestOverlap( overlapInfo &overlap, const std::string &sequence1, const std::string &sequence2 ) const;
  bool IsStaggered( const overlapInfo& overlap ) const;
  void PrintOverlap( const Sequence &seq1, const Sequence &seq2, const overlapInfo &overlap ) const;
};
