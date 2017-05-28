#pragma once

#include <string>
#include <cassert>
#include <ctype.h>
#include <numeric>
#include <map>

static bool IsBlank( const std::string &str ) {
  return str.empty() || std::all_of( str.begin(), str.end(), isspace );
}

static void UpcaseString( std::string &str ) {
  for( auto &ch : str )
    if( ch >= 97 && ch <= 122 )
      ch &= ~0x20;
}

static const int NUC_TOTAL = 16; // effectively
static const int NUC_MATRIX_SIZE = 26; // 'A'...'Z', waste some space for faster lookup
static const char NUC_MIN_ASCII = 'A';

static const int8_t NUC_MATRIX[ NUC_MATRIX_SIZE ][ NUC_MATRIX_SIZE ] = {
  {  1, -1, -1,  1,  0,  0, -1,  1,  0,  0, -1,  0,  1,  1,  0,  0,  0,  1, -1, -1, -1,  1,  1,  0, -1,  0 },
  { -1,  1,  1,  1,  0,  0,  1,  1,  0,  0,  1,  0,  1,  1,  0,  0,  0,  1,  1,  1,  1,  1,  1,  0,  1,  0 },
  { -1,  1,  1, -1,  0,  0, -1,  1,  0,  0, -1,  0,  1,  1,  0,  0,  0, -1,  1, -1, -1,  1, -1,  0,  1,  0 },
  {  1,  1, -1,  1,  0,  0,  1,  1,  0,  0,  1,  0,  1,  1,  0,  0,  0,  1,  1,  1,  1,  1,  1,  0,  1,  0 },
  {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  { -1,  1, -1,  1,  0,  0,  1, -1,  0,  0,  1,  0, -1,  1,  0,  0,  0,  1,  1, -1, -1,  1, -1,  0, -1,  0 },
  {  1,  1,  1,  1,  0,  0, -1,  1,  0,  0,  1,  0,  1,  1,  0,  0,  0,  1,  1,  1,  1,  1,  1,  0,  1,  0 },
  {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  { -1,  1, -1,  1,  0,  0,  1,  1,  0,  0,  1,  0, -1,  1,  0,  0,  0,  1,  1,  1,  1,  1,  1,  0,  1,  0 },
  {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  {  1,  1,  1,  1,  0,  0, -1,  1,  0,  0, -1,  0,  1,  1,  0,  0,  0,  1,  1, -1, -1,  1,  1,  0,  1,  0 },
  {  1,  1,  1,  1,  0,  0,  1,  1,  0,  0,  1,  0,  1,  1,  0,  0,  0,  1,  1,  1,  1,  1,  1,  0,  1,  0 },
  {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  {  1,  1, -1,  1,  0,  0,  1,  1,  0,  0,  1,  0,  1,  1,  0,  0,  0,  1,  1, -1, -1,  1,  1,  0, -1,  0 },
  { -1,  1,  1,  1,  0,  0,  1,  1,  0,  0,  1,  0,  1,  1,  0,  0,  0,  1,  1, -1, -1,  1, -1,  0,  1,  0 },
  { -1,  1, -1,  1,  0,  0, -1,  1,  0,  0,  1,  0, -1,  1,  0,  0,  0, -1, -1,  1,  1, -1,  1,  0,  1,  0 },
  { -1,  1, -1,  1,  0,  0, -1,  1,  0,  0,  1,  0, -1,  1,  0,  0,  0, -1, -1,  1,  1, -1,  1,  0,  1,  0 },
  {  1,  1,  1,  1,  0,  0,  1,  1,  0,  0,  1,  0,  1,  1,  0,  0,  0,  1,  1, -1, -1,  1,  1,  0,  1,  0 },
  {  1,  1, -1,  1,  0,  0, -1,  1,  0,  0,  1,  0,  1,  1,  0,  0,  0,  1, -1,  1,  1,  1,  1,  0,  1,  0 },
  {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
  { -1,  1,  1,  1,  0,  0, -1,  1,  0,  0,  1,  0,  1,  1,  0,  0,  0, -1,  1,  1,  1,  1,  1,  0,  1,  0 },
  {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }
};

static inline bool DoNucleotidesMatch( char nucA, char nucB ) {
  int a = nucA - NUC_MIN_ASCII;
  int b = nucB - NUC_MIN_ASCII;
  assert( a >= 0 && a < NUC_MATRIX_SIZE );
  assert( b >= 0 && b < NUC_MATRIX_SIZE );
  return NUC_MATRIX[ a ][ b ] > 0;
}

static inline char NucleotideComplement( char nuc ) {
  switch( nuc ) {
    case 'A': return 'T';
    case 'G': return 'C';
    case 'C': return 'G';
    case 'T': return 'A';
    case 'U': return 'A';

    case 'Y': return 'R';
    case 'R': return 'Y';
    case 'W': return 'W';
    case 'S': return 'S';
    case 'K': return 'M';
    case 'M': return 'K';

    case 'D': return 'H';
    case 'V': return 'B';
    case 'H': return 'D';
    case 'B': return 'V';
    case 'N': return 'N';
  }

  return nuc;
}

class Coverage {
public:
  Coverage( size_t totalSize )
    : mTotalSize( totalSize )
  {

  }

  size_t NumNonOverlaps() const {
    return mRanges.size();
  }

  size_t TotalSize() const {
    return mTotalSize;
  }

  size_t CoveredSize() const {
    size_t coveredSize = 0;
    for( auto &range : mRanges ) {
      coveredSize += ( range.second - range.first ) + 1;
    }
    return coveredSize;
  }

  float CoveredFraction() const {
    return float( CoveredSize() ) / float( TotalSize() );
  }

  void Add( size_t rangeMin, size_t rangeMax ) {
    assert( rangeMin <= rangeMax );

    std::map< size_t, size_t >::iterator insertIt, afterIt = mRanges.upper_bound( rangeMin );

    if( afterIt == mRanges.begin() || std::prev( afterIt )->second < rangeMin ) {
      insertIt = mRanges.insert( afterIt, std::pair< size_t, size_t >( rangeMin, rangeMax ) );
    } else {
      insertIt = std::prev( afterIt );
      if( insertIt->second < rangeMax ) {
        // Extend overlap
        insertIt->second = rangeMax;
      }
    }

    // Merge following ranges if covered fully by this one
    while( afterIt != mRanges.end() && rangeMax >= afterIt->first ) {
      insertIt->second = std::max( afterIt->second, insertIt->second );
      afterIt = mRanges.erase( afterIt );
    }
  }

private:
  size_t mTotalSize;
  std::map< size_t, size_t > mRanges; // Non overlapping ranges
};

