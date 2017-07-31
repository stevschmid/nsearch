#pragma once

#include "../Sequence.h"

namespace PairedEnd {
// TODO: Find good defaults! Seem to stringent
static const int    MERGER_DEFAULT_MIN_OVERLAP  = 16; // bases
static const double MERGER_DEFAULT_MIN_IDENTITY = 0.9;

class Merger {
public:
  Merger( int    minOverlap  = MERGER_DEFAULT_MIN_OVERLAP,
          double minIdentity = MERGER_DEFAULT_MIN_IDENTITY );
  bool Merge( Sequence& merged, const Sequence& fwd,
              const Sequence& rev ) const;

private:
  int    mMinOverlap;
  double mMinIdentity;

  typedef struct {
    size_t length;
    size_t pos1;
    size_t pos2;
  } overlapInfo;

  double ComputeOverlapScore( const char* sequence1, const char* sequence2,
                              const char* quality1, const char* quality2,
                              size_t len ) const;
  bool FindBestOverlap( overlapInfo& overlap, const Sequence& sequence1,
                        const Sequence& sequence2 ) const;
  bool IsStaggered( const overlapInfo& overlap ) const;
  void PrintOverlap( const Sequence& seq1, const Sequence& seq2,
                     const overlapInfo& overlap ) const;
};
} // namespace PairedEnd
