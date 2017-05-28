#pragma once

#include <unordered_map>
#include <set>

#include "Sequence.h"
#include "Utils.h"

struct SequenceHasher {
  size_t operator()( const std::string& sequence ) const {
    // Word has to fit
    assert( sequence.size() * 2 <= sizeof( size_t ) * 8 );

    size_t key = 0;

    for( int k = 0; k < sequence.size(); k++ ) {
      int val = 0;
      switch( sequence[ k ] ) {
        case 'A': val = 0b00; break;
        case 'C': val = 0b01; break;
        case 'T': val = 0b10; break;
        case 'G': val = 0b11; break;
        default: assert( false ); break;
      }

      key |= ( val << k * 2 );
    }

    return key;
  }
};

typedef struct {
  int score;
  int mismatches;
  double identity;
  Cigar cigar;
  std::string quality;
} Hit;

class Database {
  typedef std::pair< size_t, const Sequence* > SequenceInfo;

public:
  Database( int wordSize )
    : mWordSize( wordSize )
  {
    // Sequences have to fit in the calculated hash
    assert( wordLength * 2 <= sizeof( size_t ) * 8 );
  }

  void AddSequence( const Sequence &seq ) {
    // Save
    mSequences.push_back( seq );
    const Sequence *ref = &mSequences.back();

    // Map words
    ForEveryWordInSequence( *ref, [&]( size_t pos, std::string& word ) {
      this->mWords[ word ].push_back( SequenceInfo( pos, ref ) );
    });
  }

  SequenceList Query( const Sequence &query, int maxHits = 10 ) const {
    std::unordered_map< const Sequence*, Coverage > candidates;

    ForEveryWordInSequence( query, [&]( size_t pos, std::string &word ) {
      for( auto &seqInfo : mWords.at( word ) ) {
        size_t candidatePos = seqInfo.first;
        const Sequence *candidateSeq = seqInfo.second;

        if( candidates.find( candidateSeq ) == candidates.end() ) {
          candidates[ candidateSeq ] = Coverage( query.Length() );
        }

        Coverage &coverage = candidates[ candidateSeq ];
        coverage.Add( pos, pos + mWordSize - 1 );
      }
    });

    std::set< std::pair< Coverage, const Sequence* > > sortedCandidates;
    for( auto &c : candidates )
      sortedCandidates.insert( std::pair< Coverage, const Sequence* >( c.second, c.first ) );

    SequenceList list;
    for( auto it = sortedCandidates.rbegin(); it != sortedCandidates.rend(); ++it ) {
      /* std::cout << "Candidate " << (*it).first.CoveredFraction() << std::endl; */
      list.push_back( *(*it).second );

      if( list.size() >= maxHits )
        break;
    }
    /* std::cout << "====== " << std::endl; */

    return list;
  }

  size_t NumWords() const {
    return mWords.size();
  }

private:
  void ForEveryWordInSequence( const Sequence &seq, const std::function< void( size_t, std::string& ) > &block ) const {
    int wordLength = std::min< int >( seq.Length(), mWordSize );
    for( int i = 0; i <= seq.Length() - wordLength; i++ ) {
      std::string word = seq.sequence.substr( i, wordLength );
      ForEveryWordVariation( block, word, i );
    }
  }

  void ForEveryWordVariation( const std::function< void( size_t, std::string& ) > &block,
      std::string &word, size_t global, size_t start = 0 ) const
  {
    static const char BASES[] = "ATCG";

    // Find first ambiguous nucleotide
    size_t pos = word.find_first_not_of( BASES, start );
    if( pos != std::string::npos ) {
      bool skip = false;

      if( !skip ) {
        // Iterate through all possible bases
        char nuc = word[ pos ];
        for( int i = 0; i < sizeof( BASES ) - 1; i++ ) {
          if( DoNucleotidesMatch( nuc, BASES[ i ] ) )  {
            word[ pos ] = BASES[ i ];
            ForEveryWordVariation( block, word, global, pos + 1 );
          }
        }
      }
    } else {
      block( global, word );
    }
  }

  int mWordSize;

  SequenceList mSequences;
  std::unordered_map< std::string, std::deque< SequenceInfo >, SequenceHasher > mWords;
};
