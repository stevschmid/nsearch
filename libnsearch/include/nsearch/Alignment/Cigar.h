#pragma once

enum class CigarOp : char {
  UNKNOWN   = ' ',
  MATCH     = 'M',
  MISMATCH  = 'X',
  DELETION  = 'D',
  INSERTION = 'I',
};
using CigarOps = std::vector< CigarOp >;

class CigarEntry {
public:
  int count = 0;
  CigarOp op = CigarOp::UNKNOWN;

  CigarEntry() { }
  CigarEntry( int count, CigarOp op )
    : count( count ), op( op ) { }
};

class Cigar : public std::deque< CigarEntry > {
public:
  Cigar operator+( const Cigar &other ) const {
    Cigar ce = *this;
    for( auto &c : other )
      ce.Add( c );
    return ce;
  }

  Cigar& operator+=( const Cigar &other ) {
    for( auto &c : other )
      Add( c );
    return *this;
  }

  void Clear() {
    clear();
  }

  void Reverse() {
    std::reverse( begin(), end() );
  }

  void Add( const CigarOp &op ) {
    Add( CigarEntry( 1, op ) );
  }

  void Add( const CigarEntry &entry ) {
    if( entry.count == 0)
      return;

    if( entry.op == CigarOp::UNKNOWN )
      return;

    if( empty() ) {
      push_back( entry );
    } else {
      auto &last = *rbegin();
      if( last.op == entry.op ) {
        // merge
        last.count += entry.count;
      } else {
        push_back( entry );
      }
    }
  }

  std::string ToString() const {
    std::stringstream ss;
    for( auto &c : *this )  {
      ss << c.count << ( char )c.op;
    }
    return ss.str();
  }
};

static std::ostream &operator<<( std::ostream &os, const Cigar &cigar ) {
  return ( os << cigar.ToString() );
}
