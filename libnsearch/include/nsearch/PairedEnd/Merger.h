#pragma once

#include "../Sequence.h"

namespace PairedEnd {
// TODO: Find good defaults! Seem to stringent
static const int    MERGER_DEFAULT_MIN_OVERLAP  = 16; // bases
static const double MERGER_DEFAULT_MIN_IDENTITY = 0.9;

template< typename Alphabet >
class Merger {
public:
  Merger( const int    minOverlap  = MERGER_DEFAULT_MIN_OVERLAP,
          const double minIdentity = MERGER_DEFAULT_MIN_IDENTITY );
  bool Merge( const Sequence< Alphabet >& fwd, const Sequence< Alphabet >& rev,
              Sequence< Alphabet >* merged ) const;

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
  bool   FindBestOverlap( const Sequence< Alphabet >& sequence1,
                          const Sequence< Alphabet >& sequence2,
                          OverlapInfo*                overlap ) const;
  bool IsStaggered( const OverlapInfo& overlap ) const;
  void   PrintOverlap( const Sequence< Alphabet >& seq1,
                       const Sequence< Alphabet >& seq2,
                       const OverlapInfo&          overlap ) const;
};

} // namespace PairedEnd

#include "../../../src/PairedEnd/Merger.cpp"

