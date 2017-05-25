#include "nsearch/LocalAligner.h"

#include <ksw.h>

#include <iostream>

LocalAligner::LocalAligner( int matchScore, int mismatchScore, int gapOpenPenalty, int gapExtendPenalty )
  : mGapOpenPenalty( gapOpenPenalty ), mGapExtendPenalty( gapExtendPenalty )
{
  for( int y = 0; y < NUC_MATRIX_SIZE; y++ ) {
    for( int x = 0; x < NUC_MATRIX_SIZE; x++ ) {
      mScoringMatrix[ y * NUC_MATRIX_SIZE + x ] =
        NUC_MATRIX[ y ][ x ] > 0 ? matchScore : mismatchScore;
    }
  }
}


void LocalAligner::Align( const Sequence &query, const Sequence& target ) {
  std::string binQuery = query.sequence;
  for( auto &ch : binQuery )
    ch -= NUC_MIN_ASCII;

  std::string binTarget = target.sequence;
  for( auto &ch : binTarget )
    ch -= NUC_MIN_ASCII;

  // If one query is aligned against multiple
  // target sequences, *qry should be set to NULL during the first call and
  // freed after the last call.
  kswq_t *qry = NULL;

  kswr_t result = ksw_align(
      binQuery.length(),
      (uint8_t*)&binQuery[ 0 ],

      binTarget.length(),
      (uint8_t*)&binTarget[ 0 ],

      NUC_MATRIX_SIZE,
      mScoringMatrix,

      mGapOpenPenalty,
      mGapExtendPenalty,
      KSW_XSTART,
      &qry );

  std::cout << "score " << result.score << std::endl;
  std::cout << "score2 " << result.score2 << std::endl;
  std::cout << "query " << result.qe << " - " <<  result.qb << std::endl;
  std::cout << "target " << result.te << " - " <<  result.tb << std::endl;

  // Make sure to free stuff
  free( qry );
}
