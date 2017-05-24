#pragma once

#include <deque>
#include <string>

/**
 *
 * DNA Sequence. Allowed nucleotides: http://www.bioinformatics.org/sms/iupac.html
 *
 */
class Sequence {
public:
  Sequence();
  Sequence( const std::string &sequence );
  Sequence( const std::string &identifier, const std::string &sequence );
  Sequence( const std::string &identifier, const std::string &sequence, const std::string &quality );

  size_t Length() const;

  Sequence Subsequence( size_t pos, size_t len = std::string::npos ) const;

  Sequence operator+( const Sequence& other ) const;
  char& operator[]( size_t index );

  Sequence Reverse() const;
  Sequence Complement() const; // complements only ATCG
  Sequence ReverseComplement() const;

  std::string sequence;
  std::string identifier;
  std::string quality;
};

typedef std::deque< Sequence > SequenceList;
