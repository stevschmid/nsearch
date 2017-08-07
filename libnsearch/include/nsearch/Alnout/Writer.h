#pragma once

#include "../Database/GlobalSearch.h"
#include "../Sequence.h"

#include "../Alphabet/DNA.h"
#include "../Alphabet/Protein.h"

#include <fstream>

#define MAX_ALIGNMENT_STRING_LENGTH_LINE 60

namespace Alnout {

template <typename Alphabet>
class Writer {
private:
  std::ofstream mFile;
  std::ostream& mOutput;

  static inline char MatchSymbol( const char A, const char B );
  static inline std::string Unit();

public:
  Writer( std::ostream& output ) : mOutput( output ) {}

  Writer( const std::string& pathToFile )
      : mFile( pathToFile ), mOutput( mFile ) {}

  void operator<<(
    const std::pair< Sequence< Alphabet >, HitList< Alphabet > >&
      queryWithHits ) {
    const auto& query = queryWithHits.first;
    const auto& hits  = queryWithHits.second;

    // Output with fixed precision (sticky)
    mOutput << std::setiosflags( std::ios::fixed );

    mOutput << "Query >" << query.identifier << std::endl;
    mOutput << " %Id   TLen  Target" << std::endl;
    for( auto& hit : hits ) {
      mOutput << std::setprecision( 0 ) << std::setw( 3 )
              << ( hit.alignment.Identity() * 100.0 ) << '%' << std::setw( 7 )
              << hit.target.Length() << "  " << hit.target.identifier
              << std::endl;
    }
    mOutput << std::endl;

    for( auto& hit : hits ) {
      auto queryLen  = std::to_string( query.Length() );
      auto targetLen = std::to_string( hit.target.Length() );
      auto maxLen    = std::max( queryLen.size(), targetLen.size() );

      mOutput << " Query" << std::setw( maxLen + 1 )
              << std::to_string( query.Length() ) << Unit()
              << " >" << query.identifier << std::endl;
      mOutput << "Target" << std::setw( maxLen + 1 )
              << std::to_string( hit.target.Length() ) << Unit()
              << " >" << hit.target.identifier << std::endl;

      size_t numCols, numMatches, numGaps;
      auto   lines = ExtractAlignmentLines( query, hit.target, hit.alignment,
                                          &numCols, &numMatches, &numGaps );

      mOutput << std::endl;
      for( auto& line : lines ) {
        auto padLen = std::max( std::to_string( lines.back().qe ).size(),
                                std::to_string( lines.back().te ).size() );

        mOutput << "Qry " << std::setw( padLen ) << line.qs
                << " "
                << line.q << " " << line.qe << std::endl;

        mOutput << std::string( 5 + padLen, ' ' ) << line.a << std::endl;

        mOutput << "Tgt " << std::setw( padLen ) << line.ts
                << " "
                << line.t << " " << line.te << std::endl;

        mOutput << std::endl;
      }

      float identity  = float( numMatches ) / float( numCols );
      float gapsRatio = float( numGaps ) / float( numCols );
      mOutput << numCols << " cols, " << numMatches << " ids ("
              << std::setprecision( 1 ) << ( 100.0f * identity ) << "%), "
              << numGaps << " gaps (" << std::setprecision( 1 )
              << ( 100.0f * gapsRatio ) << "%)" << std::endl;

      mOutput << std::endl;
    }
  }

private:
  using AlignmentLine = struct {
    size_t      qs, qe;
    std::string q;

    size_t      ts, te;
    std::string t;

    std::string a;
  };
  using AlignmentLines = std::deque< AlignmentLine >;

  static AlignmentLines ExtractAlignmentLines( const Sequence< Alphabet >& query,
                                               const Sequence< Alphabet >& target,
                                               const Cigar&    alignment,
                                               size_t* outNumCols    = NULL,
                                               size_t* outNumMatches = NULL,
                                               size_t* outNumGaps    = NULL ) {
    Cigar cigar = alignment;

    size_t queryStart  = 0;
    size_t targetStart = 0;

    // Dont take left terminal gap into account
    if( !cigar.empty() ) {
      const auto& fce = cigar.front();
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
      const auto& bce = cigar.back();
      if( bce.op == CigarOp::DELETION ) {
        cigar.pop_back();
      } else if( bce.op == CigarOp::INSERTION ) {
        cigar.pop_back();
      }
    }

    bool   match;
    size_t numMatches = 0;
    size_t numCols    = 0;
    size_t numGaps    = 0;

    size_t qcount = queryStart;
    size_t tcount = targetStart;

    AlignmentLine line;
    line.qs = queryStart + 1;
    line.ts = targetStart + 1;

    AlignmentLines lines;

    for( auto& c : cigar ) {
      for( int i = 0; i < c.count; i++ ) {
        switch( c.op ) {
          case CigarOp::INSERTION:
            line.t += '-';
            line.q += query[ qcount++ ];
            line.a += ' ';
            numGaps++;
            break;

          case CigarOp::DELETION:
            line.q += '-';
            line.t += target[ tcount++ ];
            line.a += ' ';
            numGaps++;
            break;

          case CigarOp::MATCH: // specific for match
            numMatches++;
          case CigarOp::MISMATCH: // match and mismatch
            line.q += query[ qcount++ ];
            line.t += target[ tcount++ ];
            {
              const char a = line.q.back(), b = line.t.back();
              line.a += MatchSymbol( a, b );
            }
            break;

          default:
            break;
        }

        numCols++;
        if( numCols % MAX_ALIGNMENT_STRING_LENGTH_LINE == 0 ) {
          line.qe = qcount;
          line.te = tcount;
          lines.push_back( line );

          // Start new line
          line    = AlignmentLine();
          line.qs = qcount + 1;
          line.ts = tcount + 1;
        }
      }
    }

    if( !line.a.empty() ) {
      line.qe = qcount;
      line.te = tcount;
      lines.push_back( line );
    }

    if( outNumCols )
      *outNumCols = numCols;

    if( outNumMatches )
      *outNumMatches = numMatches;

    if( outNumGaps )
      *outNumGaps = numGaps;

    return lines;
  }

}; // Writer

template<>
inline char Writer< DNA >::MatchSymbol( const char A, const char B ) {
  return MatchPolicy< DNA >::Match( A, B ) ? '|' : ' ';
}

template<>
inline std::string Writer< DNA >::Unit() {
  return "nt";
}

template<>
inline char Writer< Protein >::MatchSymbol( const char A, const char B ) {
  auto score = ScorePolicy< Protein >::Score( A, B );
  if( score >= 4 )
    return '|';
  if( score >= 2 )
    return ':';
  if( score > 0 )
    return '.';
  return ' ';
}

template<>
inline std::string Writer< Protein >::Unit() {
  return "aa";
}

} // namespace Alnout
