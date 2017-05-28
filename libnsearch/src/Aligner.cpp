#include "nsearch/Aligner.h"

#include <ksw.h>

std::string BinarifySequence( const Sequence &seq ) {
  std::string binSeq = seq.sequence;
  for( auto &ch : binSeq )
    ch -= NUC_MIN_ASCII;
  return binSeq;
}


Aligner::Aligner( int matchScore, int mismatchScore, int gapOpenPenalty, int gapExtendPenalty )
  : mGapOpenPenalty( gapOpenPenalty ), mGapExtendPenalty( gapExtendPenalty )
{
  for( int y = 0; y < NUC_MATRIX_SIZE; y++ ) {
    for( int x = 0; x < NUC_MATRIX_SIZE; x++ ) {
      mScoringMatrix[ y * NUC_MATRIX_SIZE + x ] =
        NUC_MATRIX[ y ][ x ] > 0 ? matchScore : mismatchScore;
    }
  }
}

int Aligner::LocalAlign( const Sequence &query, const Sequence& target, Alignment *alignment, QueryProfileCache *queryProfile ) const {
  std::string binQuery = BinarifySequence( query );
  std::string binTarget = BinarifySequence( target );

  // If one query is aligned against multiple
  // target sequences, *qry should be set to NULL during the first call and
  // freed after the last call.
  kswq_t *qry = NULL;
  if( queryProfile ) {
    qry = ( kswq_t* )(*queryProfile).get();
  }

  int extraFlags = 0;
  if( alignment )  {
    extraFlags |= KSW_XSTART;
  }

  kswr_t result = ksw_align(
      binQuery.length(),
      (uint8_t*)&binQuery[ 0 ],

      binTarget.length(),
      (uint8_t*)&binTarget[ 0 ],

      NUC_MATRIX_SIZE,
      mScoringMatrix,

      mGapOpenPenalty,
      mGapExtendPenalty,
      extraFlags,
      &qry );

  if( queryProfile ) {
    if( (*queryProfile).get() != qry ) {
      *queryProfile = QueryProfileCache( qry, [=]( void *ptr ) { free( ptr ); } );
    }
  } else {
    free( qry );
  }

  if( alignment ) {
    int targetStart = result.tb;
    int targetLength = ( result.te - result.tb ) + 1;
    int queryStart = result.qb;
    int queryLength = ( result.qe - result.qb ) + 1;

    int score = GlobalAlign( query.Subsequence( queryStart, queryLength ),
        target.Subsequence( targetStart, targetLength ),
        alignment );

    int gapHead = targetStart;
    int gapTail = target.Length() - ( targetStart + targetLength );

    if( gapHead > 0 ) alignment->cigar.push_front( CigarPair( gapHead, CIGAR_SOFT_CLIP ) );
    if( gapTail > 0 ) alignment->cigar.push_back( CigarPair( gapTail, CIGAR_SOFT_CLIP ) );
  }

  return result.score;
}

int Aligner::GlobalAlign( const Sequence &query, const Sequence& target, Alignment *alignment ) const {
  std::string binQuery = BinarifySequence( query );
  std::string binTarget = BinarifySequence( target );

  int queryLen = binQuery.length();
  int targetLen = binTarget.length();
  int bandWidth = std::max< int >( 50, 1.5 * abs( queryLen - targetLen ) );

  int numCigar;
  uint32_t *cigar;

  // Toggle backtracking depending if we want Alignment or just score
  int *numCigarArg = NULL;
  uint32_t **cigarArg = NULL;
  if( alignment ) {
    numCigarArg = &numCigar;
    cigarArg = &cigar;
  }

  int score = ksw_global(
      queryLen,
      (uint8_t*)&binQuery[ 0 ],

      targetLen,
      (uint8_t*)&binTarget[ 0 ],

      NUC_MATRIX_SIZE,
      mScoringMatrix,

      mGapOpenPenalty,
      mGapExtendPenalty,

      bandWidth,
      numCigarArg,
      cigarArg );

  if( alignment ) {
    alignment->cigar.clear();
    for( int i = 0; i < numCigar; i++ ) {
      CigarPair cp;

      int num = cigar[ i ] >> 4;
      cp.first = num;

      int op = cigar[ i ] & 0b1111;
      switch( op ) {
        case 0: cp.second = CIGAR_MATCH; break;
        case 1: cp.second = CIGAR_INSERTION; break;
        case 2: cp.second = CIGAR_DELETION; break;
      }

      alignment->cigar.push_back( std::move( cp ) );
    }

    // delete cigar
    free( cigar );
  }

  return score;
}

void PrintAlignment( const Alignment &aln, const Sequence &query, const Sequence &target, std::ostream &os ) {
  int qcount = 0;
  int tcount = 0;

  std::string q;
  std::string t;

  for( auto &p : aln.cigar ) {
    int len = p.first;
    CigarOperation op = p.second;

    for( int i = 0; i < len; i++ ) {
      switch( op ) {
        case 'I':
          t += '.';
          q += query[ qcount++ ];
          break;

        case 'S':
        case 'H':
        case 'D':
          q += '.';
          t += target[ tcount++ ];
          break;

        case 'M':
          q += query[ qcount++ ];
          t += target[ tcount++ ];
          break;

        default:
          break;
      }
    }
  }

  os << "REF " << t << std::endl;
  os << "QRY " << q << std::endl;
}
