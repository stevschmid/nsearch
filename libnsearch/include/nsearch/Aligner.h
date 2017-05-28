#pragma once

#include "Sequence.h"
#include "Utils.h"

#include <iostream>
#include <sstream>

typedef enum {
  CIGAR_MATCH = 'M',
  CIGAR_DELETION = 'D',
  CIGAR_INSERTION = 'I',
  CIGAR_SOFT_CLIP = 'S',
  CIGAR_HARD_CLIP = 'H',
} CigarOperation;

typedef struct std::pair< int, CigarOperation > CigarPair;
typedef std::deque< CigarPair > Cigar;

typedef struct {
  int score;
  Cigar cigar;

  std::string cigarString() {
    std::stringstream str;
    for( auto &p : cigar ) {
      str << (int)p.first << (char)p.second;
    }
    return str.str();
  }
} Alignment;

typedef std::shared_ptr< void > QueryProfileCache;

// Alignment wrapper for KSW
class Aligner {
public:
  // Gap of length L costs -(gapOpenPenalty + L * gapExtendPenalty)
  Aligner( int matchScore = 1, int mismatchScore = -2,
      int gapOpenPenalty = 10, int gapExtendPenalty = 1 );

  // Cache QueryProfileCache for same query to speed up local alignments
  int LocalAlign( const Sequence &query, const Sequence& target, Alignment *klignment = NULL, QueryProfileCache *queryProfile = NULL ) const;

  int GlobalAlign( const Sequence &query, const Sequence& target, Alignment *alignment = NULL ) const;

private:
  int mGapOpenPenalty, mGapExtendPenalty;
  int8_t mScoringMatrix[ NUC_MATRIX_SIZE * NUC_MATRIX_SIZE ];
};

extern void PrintAlignment( const Alignment &aln, const Sequence &query, const Sequence &target, std::ostream &os = std::cout );
