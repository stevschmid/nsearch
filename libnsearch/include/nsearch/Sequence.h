#pragma once

#include <cassert>
#include <deque>
#include <iostream>
#include <string>

/**
 *
 * DNA Sequence. Allowed nucleotides:
 * http://www.bioinformatics.org/sms/iupac.html
 *
 */
class Sequence {
public:
  Sequence();
  Sequence( const Sequence& sequence );
  Sequence( Sequence&& sequence );
  Sequence& operator=( const Sequence& other );
  Sequence( const std::string& sequence );
  Sequence( const char* sequence );
  Sequence( const std::string& identifier, const std::string& sequence );
  Sequence( const std::string& identifier, const std::string& sequence,
            const std::string& quality );

  size_t Length() const;

  Sequence Subsequence( const size_t pos,
                        const size_t len = std::string::npos ) const;

  Sequence operator+( const Sequence& other ) const;
  char& operator[]( const size_t index );
  char operator[]( const size_t index ) const;
  bool operator==( const Sequence& other ) const;
  bool operator!=( const Sequence& other ) const;

  Sequence Reverse() const;
  Sequence Complement() const;
  Sequence ReverseComplement() const;

  std::string sequence;
  std::string identifier;
  std::string quality;
};

static std::ostream& operator<<( std::ostream& os, const Sequence& seq ) {
  if( !seq.identifier.empty() )
    os << ">" << seq.identifier << std::endl;
  if( !seq.sequence.empty() )
    os << " " << seq.sequence << std::endl;
  if( !seq.quality.empty() )
    os << " " << seq.quality << std::endl;
  return os;
}

typedef std::deque< Sequence > SequenceList;

// Hash function (only for unambiguous DNA sequences)
namespace std
{
  template <>
  struct hash< Sequence > {
    size_t operator()( const Sequence &seq ) const {
      size_t key = 0;

      for( int k = 0; k < seq.sequence.size(); k++ ) {
        int val = 0;
        switch( seq.sequence[ k ] ) {
          case 'A': val = 0b00; break;
          case 'C': val = 0b01; break;
          case 'U':
          case 'T': val = 0b10; break;
          case 'G': val = 0b11; break;
          default: assert( false ); break;
        }

        key |= ( val << ( ( k * 2 ) % ( sizeof( size_t ) * 8 ) ) );
      }

      return key;
    }
  };
}
