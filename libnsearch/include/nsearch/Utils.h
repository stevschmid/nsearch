#pragma once

#include <string>
#include <cassert>
#include <ctype.h>

static std::string trim( const std::string &tstr ) {
  std::string str = tstr;
  str.erase( str.begin(), std::find_if_not(str.begin(), str.end(), [](char c){ return std::isspace(c); }) );
  str.erase( std::find_if_not(str.rbegin(), str.rend(), [](char c){ return std::isspace(c); }).base(), str.end() );
  return str;
}

static std::string toupper( const std::string &ustr ) {
  std::string str = ustr;
  for( auto &ch : str )
    ch = toupper( ch );
  return str;
}

static const int NUC_MATRIX_SIZE = 26; // 'A'...'Z'
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

static inline bool AreNucleotidesMatching( char nucA, char nucB ) {
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
