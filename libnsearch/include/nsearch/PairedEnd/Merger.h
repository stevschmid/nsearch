#pragma once

#include "../Sequence.h"

namespace PairedEnd {
// TODO: Find good defaults! Seem to stringent
static const int    MERGER_DEFAULT_MIN_OVERLAP  = 16; // bases
static const double MERGER_DEFAULT_MIN_IDENTITY = 0.9;

class Merger {
public:
  Merger( const int    minOverlap  = MERGER_DEFAULT_MIN_OVERLAP,
          const double minIdentity = MERGER_DEFAULT_MIN_IDENTITY );
  bool Merge( const Sequence& fwd, const Sequence& rev,
              Sequence* merged ) const;

private:
  int    mMinOverlap;
  double mMinIdentity;

  typedef struct {
    size_t length;
    size_t pos1;
    size_t pos2;
  } OverlapInfo;

  double ComputeOverlapScore( const char* sequence1, const char* sequence2,
                              const char* quality1, const char* quality2,
                              const size_t len ) const;
  bool FindBestOverlap( const Sequence& sequence1, const Sequence& sequence2,
                        OverlapInfo* overlap ) const;
  bool IsStaggered( const OverlapInfo& overlap ) const;
  void PrintOverlap( const Sequence& seq1, const Sequence& seq2,
                     const OverlapInfo& overlap ) const;
};
} // namespace PairedEnd
