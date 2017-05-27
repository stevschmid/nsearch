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

typedef struct std::pair< CigarOperation, int > CigarPair;
typedef std::deque< CigarPair > Cigar;

typedef struct {
  int score;
  Cigar cigar;

  std::string cigarString() {
    std::stringstream str;
    for( auto &p : cigar ) {
      str << (int)p.second << (char)p.first;
    }
    return str.str();
  }
} Alignment;

typedef struct {
  size_t targetStart, targetLength;
  size_t queryStart, queryLength;
} LocalAlignmentInfo;

typedef std::shared_ptr< void > QueryProfileCache;

// Alignment wrapper for KSW
class Aligner {
public:
  // Gap of length L costs -(gapOpenPenalty + L * gapExtendPenalty)
  Aligner( int matchScore = 1, int mismatchScore = -2,
      int gapOpenPenalty = 10, int gapExtendPenalty = 1 );

  // Local alignment is blazing fast but we don't get detailed alignment (no cigar)
  // Cache QueryProfileCache for same query to speed up local alignments
  int LocalAlign( const Sequence &query, const Sequence& target, LocalAlignmentInfo *info = NULL, QueryProfileCache *queryProfile = NULL ) const;
  int ComputeLocalAlignment( Alignment &alignment, const Sequence &query, const Sequence& target, const LocalAlignmentInfo &info ) const;

  // Global alignment is slower but we get detailed alignment
  int GlobalAlign( const Sequence &query, const Sequence& target, Alignment *alignment = NULL ) const;

private:
  int mGapOpenPenalty, mGapExtendPenalty;
  int8_t mScoringMatrix[ NUC_MATRIX_SIZE * NUC_MATRIX_SIZE ];
};

extern void PrintAlignment( const Alignment &aln, const Sequence &query, const Sequence &target, std::ostream &os = std::cout );
