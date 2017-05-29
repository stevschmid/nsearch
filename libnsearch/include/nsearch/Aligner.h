#pragma once

#include "Sequence.h"
#include "Utils.h"

#include <iostream>
#include <sstream>

typedef enum {
  CIGAR_MATCH     = 'M',
  CIGAR_DELETION  = 'D',
  CIGAR_INSERTION = 'I',
} CigarOperation;

typedef struct CigarItem_s {
  int length;
  CigarOperation operation;

  CigarItem_s() { }
  CigarItem_s( int length, CigarOperation operation )
    : length( length ), operation( operation ) { }
} CigarItem;
typedef std::deque< CigarItem > Cigar;

typedef struct {
  int score = 0;
  Cigar cigar;

  int queryPos = 0;
  int targetPos = 0;

  int queryLength = 0;
  int targetLength = 0;

  std::string cigarString() {
    std::stringstream str;
    for( auto &p : cigar ) {
      str << p.length << (char)p.operation;
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
  int LocalAlign( const Sequence &query, const Sequence& target, Alignment *alignment = NULL, QueryProfileCache *queryProfile = NULL ) const;
  int GlobalAlign( const Sequence &query, const Sequence& target, Alignment *alignment = NULL ) const;

private:
  int mGapOpenPenalty, mGapExtendPenalty;
  int8_t mScoringMatrix[ NUC_MATRIX_SIZE * NUC_MATRIX_SIZE ];
};

extern void PrintAlignment( const Alignment &aln, const Sequence &query, const Sequence &target, std::ostream &os = std::cout );
