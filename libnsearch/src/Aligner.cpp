#include "nsearch/Aligner.h"

#include <iostream>
#include <sstream>

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

LocalAlignment Aligner::LocalAlign( const Sequence &query, const Sequence& target, QueryProfile *queryProfile ) const {
  std::string binQuery = BinarifySequence( query );
  std::string binTarget = BinarifySequence( target );

  // If one query is aligned against multiple
  // target sequences, *qry should be set to NULL during the first call and
  // freed after the last call.
  kswq_t *qry = NULL;
  if( queryProfile ) {
    qry = ( kswq_t* )(*queryProfile).get();
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
      KSW_XSTART,
      &qry );

  LocalAlignment aln;
  aln.score = result.score;
  aln.targetStart = result.tb;
  aln.targetLength = ( result.te - result.tb ) + 1;
  aln.queryStart = result.qb;
  aln.queryLength = ( result.qe - result.qb ) + 1;

  if( queryProfile ) {
    if( (*queryProfile).get() != qry ) {
      *queryProfile = QueryProfile( qry, [=]( void *ptr ) { free( ptr ); } );
    }
  } else {
    free( qry );
  }

  return aln;
}

GlobalAlignment Aligner::GlobalAlign( const Sequence &query, const Sequence& target ) const {
  std::string binQuery = BinarifySequence( query );
  std::string binTarget = BinarifySequence( target );

  int queryLen = binQuery.length();
  int targetLen = binTarget.length();
  int bandWidth = std::max< int >( 50, 1.5 * abs( queryLen - targetLen ) );

  int numCigar;
  uint32_t *cigar;

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
      &numCigar,
      &cigar );

  GlobalAlignment aln;
  aln.score = global;

  for( int i = 0; i < numCigar; i++ ) {
    int num = cigar[ i ] >> 4;
    int op = cigar[ i ] & 0b1111;
    CigarPair cp;
    switch( op ) {
      case 0: cp.first = CIGAR_MATCH; break;
      case 1: cp.first = CIGAR_INSERTION; break;
      case 2: cp.first = CIGAR_DELETION; break;
    }
    cp.second = num;
    aln.cigar.push_back( std::move( cp ) );
  }

  // delete cigar
  free( cigar );

  return aln;
}

std::string CigarAsString( const Cigar &cigar ) {
  std::stringstream str;
  for( auto &p : cigar ) {
    str << (int)p.second << (char)p.first;
  }
  return str.str();
}

void PrettyPrintGlobalAlignment( const Sequence &query, const Sequence &target, const GlobalAlignment &aln, std::ostream &os ) {
  int qcount = 0;
  int tcount = 0;

  std::string q;
  std::string t;

  for( auto &p : aln.cigar ) {
    CigarOperation op = p.first;
    int len = p.second;

    for( int i = 0; i < len; i++ ) {
      switch( op ) {
        case 'I':
          t += '.';
          q += query[ qcount++ ];
          break;

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
