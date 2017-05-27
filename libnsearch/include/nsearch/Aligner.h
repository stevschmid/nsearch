#pragma once

#include "Sequence.h"
#include "Utils.h"

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
} GlobalAlignment;

typedef struct {
  int score;
  size_t targetStart, targetLength;
  size_t queryStart, queryLength;
} LocalAlignment;

typedef std::shared_ptr< void > QueryProfile;

class Aligner {
public:
  // Gap of length L costs -(gapOpenPenalty + L * gapExtendPenalty)
  Aligner( int matchScore = 1, int mismatchScore = -2,
      int gapOpenPenalty = 10, int gapExtendPenalty = 1 );

  // Local alignment is blazing fast but we don't get detailed alignment (no cigar)
  // Cache QueryProfile for saem query to speed up local alignments
  LocalAlignment LocalAlign( const Sequence &query, const Sequence& target, QueryProfile *queryProfile = NULL ) const;

  // Global alignment is slower but we get detailed alignment
  GlobalAlignment GlobalAlign( const Sequence &query, const Sequence& target ) const;

private:
  int mGapOpenPenalty, mGapExtendPenalty;
  int8_t mScoringMatrix[ NUC_MATRIX_SIZE * NUC_MATRIX_SIZE ];
};

extern void PrettyPrintGlobalAlignment( const Sequence &query, const Sequence &target, const GlobalAlignment &aln, std::ostream &os = std::cout );
extern std::string CigarAsString( const Cigar &cigar );

