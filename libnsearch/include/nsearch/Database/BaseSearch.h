#pragma once

#include "nsearch/Sequence.h"
#include "nsearch/Alignment/Cigar.h"

#include <deque>
#include <vector>

class Database;

class BaseSearch {
public:
  using Hit = struct {
    Sequence target;
    Cigar alignment;
  };
  using HitList = std::deque< Hit >;

  BaseSearch( const Database &db ) : mDB( db ) { }
  virtual HitList Query( const Sequence &query ) = 0;

private:
  const Database &mDB;
};

/* static std::ostream& operator<<( std::ostream &os, const BaseSearch::Match &result ) { */
/*   std::string q; */
/*   std::string t; */
/*   std::string a; */

/*   const Sequence &query = result.query; */
/*   const Sequence &target = result.target; */

/*   Cigar cigar = result.alignment; */

/*   size_t queryStart = 0; */
/*   size_t targetStart = 0; */

/*   // Dont take left terminal gap into account */
/*   if( !cigar.empty() ) { */
/*     const CigarEntry &fce = cigar.front(); */
/*     if( fce.op == CigarOp::DELETION ) { */
/*       targetStart = fce.count; */
/*       cigar.pop_front(); */
/*     } else if( fce.op == CigarOp::INSERTION ) { */
/*       queryStart = fce.count; */
/*       cigar.pop_front(); */
/*     } */
/*   } */

/*   // Don't take right terminal gap into account */
/*   if( !cigar.empty() ) { */
/*     const CigarEntry &bce = cigar.back(); */
/*     if( bce.op == CigarOp::DELETION ) { */
/*       cigar.pop_back(); */
/*     } else if( bce.op == CigarOp::INSERTION ) { */
/*       cigar.pop_back(); */
/*     } */
/*   } */

/*   bool match; */
/*   size_t numMatches = 0; */
/*   size_t numCols = 0; */

/*   size_t qcount = queryStart; */
/*   size_t tcount = targetStart; */

/*   bool correct = true; */

/*   for( auto &c : cigar ) { */
/*     for( int i = 0; i < c.count; i++ ) { */
/*       switch( c.op ) { */
/*         case CigarOp::INSERTION: */
/*           t += '-'; */
/*           q += query[ qcount++ ]; */
/*           a += ' '; */
/*           break; */

/*         case CigarOp::DELETION: */
/*           q += '-'; */
/*           t += target[ tcount++ ]; */
/*           a += ' '; */
/*           break; */

/*         case CigarOp::MATCH: */
/*           numMatches++; */
/*           q += query[ qcount++ ]; */
/*           t += target[ tcount++ ]; */
/*           { */
/*             bool match = DoNucleotidesMatch( q.back(), t.back() ); */
/*             if( !match ) { */
/*               correct = false; */

/*               a += 'X'; */
/*             } else { */
/*               a += '|'; */
/*             } */
/*           } */
/*           break; */

/*         case CigarOp::MISMATCH: */
/*           a += ' '; */
/*           q += query[ qcount++ ]; */
/*           t += target[ tcount++ ]; */
/*           break; */

/*         default: */
/*           break; */
/*       } */

/*       numCols++; */
/*     } */
/*   } */

/*   os << "Query " << std::string( 11, ' ' ) << ">" << query.identifier << std::endl; */
/*   os << std::endl; */
/*   os << std::setw( 15 ) << queryStart + 1 << " " << q << " " << qcount << std::endl; */
/*   os << std::string( 16, ' ' ) << a << std::endl; */
/*   os << std::setw( 15 ) << targetStart + 1 << " " << t << " " << tcount << std::endl; */
/*   os << std::endl; */
/*   os << "Target " << std::string( 11, ' ' ) << ">" << target.identifier << std::endl; */

/*   os << std::endl; */
/*   float identity = float( numMatches ) / float( numCols ); */
/*   os <<  numCols << " cols, " << numMatches << " ids (" << ( 100.0f * identity ) << "%)" << std::endl; */

/*   if( !correct ) { */
/*     os << "!!!INVALID ALIGNMENT!!!" << std::endl; */
/*   } */

/*   return os; */
/* } */
