#pragma once

#include <deque>
#include <vector>
#include <iostream>
#include <sstream>

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

  std::string ToFullAlignmentString( const Sequence &query, const Sequence &target, bool *correct = NULL ) const {
    std::string q;
    std::string t;
    std::string a;

    if( correct )
      *correct = true;

    size_t queryStart = 0;
    size_t targetStart = 0;

    Cigar cigar = *this;

    // Dont take left terminal gap into account
    if( !cigar.empty() ) {
      const CigarEntry &fce = cigar.front();
      if( fce.op == CigarOp::DELETION ) {
        targetStart = fce.count;
        cigar.pop_front();
      } else if( fce.op == CigarOp::INSERTION ) {
        queryStart = fce.count;
        cigar.pop_front();
      }
    }

    // Don't take right terminal gap into account
    if( !cigar.empty() ) {
      const CigarEntry &bce = cigar.back();
      if( bce.op == CigarOp::DELETION ) {
        cigar.pop_back();
      } else if( bce.op == CigarOp::INSERTION ) {
        cigar.pop_back();
      }
    }

    bool match;
    size_t numMatches = 0;
    size_t numCols = 0;

    size_t qcount = queryStart;
    size_t tcount = targetStart;

    for( auto &c : cigar ) {
      for( int i = 0; i < c.count; i++ ) {
        switch( c.op ) {
          case CigarOp::INSERTION:
            t += '-';
            q += query[ qcount++ ];
            a += ' ';
            break;

          case CigarOp::DELETION:
            q += '-';
            t += target[ tcount++ ];
            a += ' ';
            break;

          case CigarOp::MATCH:
            numMatches++;
            q += query[ qcount++ ];
            t += target[ tcount++ ];
            {
              bool match = DoNucleotidesMatch( q.back(), t.back() );
              if( !match ) {
                if( correct )
                  *correct = false;

                a += 'X';
              } else {
                a += '|';
              }
            }
            break;

          case CigarOp::MISMATCH:
            a += ' ';
            q += query[ qcount++ ];
            t += target[ tcount++ ];
            break;

          default:
            break;
        }

        numCols++;
      }
    }

    std::stringstream ss;

    ss << "Query " << std::string( 11, ' ' ) << ">" << query.identifier << std::endl;
    ss << std::endl;
    ss << std::setw( 15 ) << queryStart + 1 << " " << q << " " << qcount << std::endl;
    ss << std::string( 16, ' ' ) << a << std::endl;
    ss << std::setw( 15 ) << targetStart + 1 << " " << t << " " << tcount << std::endl;
    ss << std::endl;
    ss << "Target " << std::string( 11, ' ' ) << ">" << target.identifier << std::endl;

    ss << std::endl;
    float identity = float( numMatches ) / float( numCols );
    ss <<  numCols << " cols, " << numMatches << " ids (" << ( 100.0f * identity ) << "%)" << std::endl;

    return ss.str();
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
