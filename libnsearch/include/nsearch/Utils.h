#pragma once

#include <string>
#include <cassert>
#include <ctype.h>
#include <numeric>

static bool IsBlank( const std::string& str ) {
  return str.empty() || std::all_of( str.begin(), str.end(), isspace );
}

static void UpcaseString( std::string& str ) {
  for( auto& ch : str )
    if( ch >= 97 && ch <= 122 ) // upcase
      ch &= ~0x20;
}

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

static inline bool DoNucleotidesMatch( const char nucA, const char nucB ) {
  int a = nucA - NUC_MIN_ASCII;
  int b = nucB - NUC_MIN_ASCII;
  assert( a >= 0 && a < NUC_MATRIX_SIZE );
  assert( b >= 0 && b < NUC_MATRIX_SIZE );
  return NUC_MATRIX[ a ][ b ] > 0;
}

static inline char NucleotideComplement( const char nuc ) {
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
