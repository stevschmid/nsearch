#include "nsearch/Utils.h"

#include <cassert>
#include <ctype.h>
#include <iostream>

template < typename A >
Sequence< A >::Sequence() : Sequence( "", "", "" ) {}

template < typename A >
Sequence< A >::Sequence( const Sequence& sequence )
    : sequence( sequence.sequence ), identifier( sequence.identifier ),
      quality( sequence.quality ) {}

template < typename A >
Sequence< A >::Sequence( Sequence< A >&& sequence )
    : sequence( std::move( sequence.sequence ) ),
      identifier( std::move( sequence.identifier ) ),
      quality( std::move( sequence.quality ) ) {}

template < typename A >
Sequence< A >& Sequence< A >::operator=( const Sequence< A >& other ) {
  sequence   = other.sequence;
  identifier = other.identifier;
  quality    = other.quality;
  return *this;
}

template < typename A >
Sequence< A >::Sequence( const std::string& sequence )
    : Sequence( "", sequence, "" ) {}

template < typename A >
Sequence< A >::Sequence( const char* sequence )
    : Sequence( "", sequence, "" ) {}

template < typename A >
Sequence< A >::Sequence(
  const std::string&                               identifier, const std::basic_string< typename A::CharType >& sequence )
    : Sequence( identifier, sequence, "" ) {}

template < typename A >
Sequence< A >::Sequence(
  const std::string&                               identifier,
  const std::basic_string< typename A::CharType >& sequence,
  const std::string&                               quality )
    : identifier( identifier ), sequence( sequence ), quality( quality ) {}

template < typename A >
size_t Sequence< A >::Length() const {
  return sequence.length();
}

template < typename A >
Sequence< A > Sequence< A >::Subsequence( const size_t pos,
                                          const size_t len_ ) const {
  size_t len = len_;
  if( len == std::string::npos ) {
    len = Length() - pos;
  }

  return Sequence( identifier,
                   pos < sequence.length() ? sequence.substr( pos, len ) : "",
                   pos < quality.length() ? quality.substr( pos, len ) : "" );
}

template < typename A >
Sequence< A > Sequence< A >::operator+( const Sequence< A >& other ) const {
  return Sequence< A >( identifier, sequence + other.sequence,
                        quality + other.quality );
}

template < typename A >
char& Sequence< A >::operator[]( const size_t index ) {
  assert( index >= 0 && index < sequence.size() );
  return sequence[ index ];
}

template < typename A >
char Sequence< A >::operator[]( const size_t index ) const {
  assert( index >= 0 && index < sequence.size() );
  return sequence[ index ];
}

template < typename A >
bool Sequence< A >::operator==( const Sequence< A >& other ) const {
  return !( *this != other );
}

template < typename A >
bool Sequence< A >::operator!=( const Sequence< A >& other ) const {
  if( Length() != other.Length() )
    return true;

  auto tit = ( *this ).sequence.begin();
  auto oit = other.sequence.begin();
  while( tit != ( *this ).sequence.end() && oit != other.sequence.end() ) {
    if( !A::Match( *tit, *oit ) )
      return true;

    ++tit;
    ++oit;
  }

  return false;
}

template < typename A >
Sequence< A > Sequence< A >::Complement() const {
  Sequence complement = *this;

  for( char& ch : complement.sequence ) {
    ch = NucleotideComplement( ch );
  }

  return complement;
}

template < typename A >
Sequence< A > Sequence< A >::Reverse() const {
  Sequence rev = *this;
  // Reverse sequence and quality
  std::reverse( rev.sequence.begin(), rev.sequence.end() );
  std::reverse( rev.quality.begin(), rev.quality.end() );
  return rev;
}

template < typename A >
Sequence< A > Sequence< A >::ReverseComplement() const {
  return Reverse().Complement();
}
