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

/*   kswr_t result = ksw_align( */
/*       binQuery.length(), */
/*       (uint8_t*)&binQuery[ 0 ], */

/*       binTarget.length(), */
/*       (uint8_t*)&binTarget[ 0 ], */

/*       NUC_MATRIX_SIZE, */
/*       mScoringMatrix, */

/*       mGapOpenPenalty, */
/*       mGapExtendPenalty, */
/*       KSW_XSTART, */
/*       &qry ); */

/*   std::cout << "score " << result.score << std::endl; */
/*   std::cout << "score2 " << result.score2 << std::endl; */
/*   std::cout << "query " << result.qe << " - " <<  result.qb << std::endl; */
/*   std::cout << "target " << result.te << " - " <<  result.tb << std::endl; */
/*   std::cout << "lenT " << result.te - result.tb + 1 << std::endl; */
/*   std::cout << "lenQ " << result.qe - result.qb + 1 << std::endl; */

/*   std::cout << std::string( result.tb, ' ') << query.sequence << std::endl; */
/*   std::cout << std::string( std::max( result.qb, result.tb ), ' ' ); */
/*   for( int i = 0; i < result.qe - result.qb + 1; i++ ) { */
/*     std::cout << '.'; */
/*     /1* if( seq1.sequence[ overlap.pos1 + i ] == seq2.sequence[ overlap.pos2 + i ] ) { *1/ */
/*     /1*   std::cout << '|'; *1/ */
/*     /1* } else { *1/ */
/*     /1*   std::cout << ' '; *1/ */
/*     /1* } *1/ */
/*   } */
/*   std::cout << std::endl; */
/*   std::cout << std::string( result.qb, ' ' ) << target.sequence << std::endl; */

  int queryLen = binQuery.length();
  int targetLen = binTarget.length();
  int bandWidth = std::max< int >( 50, 1.5 * abs( queryLen - targetLen ) );

  int numCigars;
  uint32_t *cigars;

  int global = ksw_global(
      queryLen,
      (uint8_t*)&binQuery[ 0 ],

      targetLen,
      (uint8_t*)&binTarget[ 0 ],

      NUC_MATRIX_SIZE,
      mScoringMatrix,

      mGapOpenPenalty,
      mGapExtendPenalty,

      bandWidth,
      &numCigars,
      &cigars );

  /* std::cout << global << std::endl; */
  /* std::cout << numCigars << std::endl; */

  int qcount = 0;
  int tcount = 0;

  std::string q;
  std::string t;

  for( int i = 0; i < numCigars; i++ ) {
    int len = cigars[ i ] >> 4;
    int op = cigars[ i ] & 0b1111;
    char opchar = '?';
    switch( op ) {
      case 0:
        opchar = 'M';
        break;

      case 1:
        opchar = 'I';
        break;

      case 2:
        opchar = 'D';
        break;
    }

    for( int i = 0; i < len; i++ ) {
      if( opchar == 'I' ) {
        t += '.';
        q += query[ qcount++ ];
      }
      if( opchar == 'D' ) {
        q += '.';
        t += target[ tcount++ ];
      }
      if( opchar == 'M' ) {
        q += query[ qcount++ ];
        t += target[ tcount++ ];
      }
    }

    printf("%d%c", len, opchar);
  }
  std::cout << std::endl;

  std::cout << "REF " << t << std::endl;
  std::cout << "QRY " << q << std::endl;


  // Make sure to free stuff
  free( cigars );
  free( qry );
}
