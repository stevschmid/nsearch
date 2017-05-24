#include "nsearch/Sequence.h"

#include <iostream>
#include <cassert>
#include <ctype.h>

Sequence::Sequence()
  : Sequence( "", "", "" )
{

}

Sequence::Sequence( const std::string &sequence )
  : Sequence( "", sequence, "" )
{

}
Sequence::Sequence( const std::string &identifier, const std::string &sequence )
  : Sequence( identifier, sequence, "" )
{

}

Sequence::Sequence( const std::string &identifier, const std::string &sequence, const std::string &quality )
  : identifier( identifier ), sequence( sequence ), quality( quality )
{

}

size_t Sequence::Length() const {
  return sequence.length();
}

Sequence Sequence::Subsequence( size_t pos, size_t len ) const {
  if( len == std::string::npos ) {
    len = Length() - pos;
  }

  return Sequence( identifier,
      sequence.substr( pos, len ),
      quality.substr( pos, len ) );
}

Sequence Sequence::operator+( const Sequence& other ) const {
  return Sequence( identifier,
      sequence + other.sequence,
      quality + other.quality );
}

char& Sequence::operator[]( size_t index ) {
  assert( index >= 0 && index < sequence.size() );
  return sequence[ index ];
}

Sequence Sequence::Complement() const {
  Sequence complement = *this;

  for( char &ch : complement.sequence ) {
    switch( toupper( ch ) ) {
      case 'A': ch = 'T'; break;
      case 'T': ch = 'A'; break;
      case 'C': ch = 'G'; break;
      case 'G': ch = 'C'; break;
    }
  }

  return complement;
}

Sequence Sequence::Reverse() const {
  Sequence rev = *this;
  // Reverse sequence and quality
  std::reverse( rev.sequence.begin(), rev.sequence.end() );
  std::reverse( rev.quality.begin(), rev.quality.end() );
  return rev;
}

Sequence Sequence::ReverseComplement() const {
  return Reverse().Complement();
}
