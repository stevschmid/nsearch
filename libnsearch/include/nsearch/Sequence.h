#pragma once

#include <deque>
#include <string>
#include <iostream>

/**
 *
 * DNA Sequence. Allowed nucleotides: http://www.bioinformatics.org/sms/iupac.html
 *
 */
class Sequence {
public:
  Sequence();
  Sequence( const std::string &sequence );
  Sequence( const char *sequence );
  Sequence( const std::string &identifier, const std::string &sequence );
  Sequence( const std::string &identifier, const std::string &sequence, const std::string &quality );

  size_t Length() const;

  Sequence Subsequence( size_t pos, size_t len = std::string::npos ) const;

  Sequence operator+( const Sequence &other ) const;
  char& operator[]( size_t index );
  char operator[]( size_t index ) const;
  bool operator==( const Sequence &other ) const;
  bool operator!=( const Sequence &other ) const;

  Sequence Reverse() const;
  Sequence Complement() const; // complements only ATCG
  Sequence ReverseComplement() const;

  std::string sequence;
  std::string identifier;
  std::string quality;
};

static std::ostream &operator<<( std::ostream &os, const Sequence &seq ) {
  return os << seq.identifier;
}

typedef std::deque< Sequence > SequenceList;
