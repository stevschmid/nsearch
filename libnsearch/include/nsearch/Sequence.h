#pragma once

#include <cassert>
#include <deque>
#include <iostream>
#include <string>

#include "Utils.h"

struct DNA {
  typedef char CharType;

  static bool Match( const char nucA, const char nucB ) {
    return DoNucleotidesMatch( nucA, nucB );
  }

  static char Complement( const char nuc ) {
    return NucleotideComplement( nuc );
  }
};
using RNA = DNA;

template < typename Alphabet >
class Sequence {
public:
  Sequence();
  Sequence( const Sequence< Alphabet >& sequence );
  Sequence( Sequence< Alphabet >&& sequence );
  Sequence< Alphabet >& operator=( const Sequence< Alphabet >& other );
  Sequence( const std::string& sequence );
  Sequence( const char* sequence );
  Sequence( const std::string&                                      identifier,
            const std::basic_string< typename Alphabet::CharType >& sequence );
  Sequence( const std::string&                                      identifier,
            const std::basic_string< typename Alphabet::CharType >& sequence,
            const std::string&                                      quality );

  size_t Length() const;

  Sequence< Alphabet >
  Subsequence( const size_t pos,
               const size_t len = std::string::npos ) const;

  Sequence< Alphabet > operator+( const Sequence< Alphabet >& other ) const;
  char&                operator[]( const size_t index );
  char                 operator[]( const size_t index ) const;
  bool                 operator==( const Sequence< Alphabet >& other ) const;
  bool                 operator!=( const Sequence< Alphabet >& other ) const;

  Sequence< Alphabet > Reverse() const;
  Sequence< Alphabet > Complement() const;
  Sequence< Alphabet > ReverseComplement() const;

  std::string identifier;
  std::string quality;

  std::basic_string< typename Alphabet::CharType > sequence;
};

template < typename Alphabet >
static std::ostream& operator<<( std::ostream& os, const Sequence< Alphabet >& seq ) {
  if( !seq.identifier.empty() )
    os << ">" << seq.identifier << std::endl;
  if( !seq.sequence.empty() )
    os << " " << seq.sequence << std::endl;
  if( !seq.quality.empty() )
    os << " " << seq.quality << std::endl;
  return os;
}

template< typename Alphabet >
using SequenceList = std::deque< Sequence< Alphabet > >;

#include "../../src/Sequence.cpp"
