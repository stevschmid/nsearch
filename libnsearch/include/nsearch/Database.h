#pragma once

#include <unordered_map>
#include <unordered_set>

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

class Database {
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
    ForEveryWordInSequence( *ref, [&]( std::string& word ) {
      this->mWordToSequences[ word ].push_back( ref );
    });
  }

  SequenceList Query( const Sequence &seq ) {
    std::unordered_set< const Sequence* > candidates;

    ForEveryWordInSequence( seq, [&]( std::string &word ) {
      for( const Sequence *candidateSeq : mWordToSequences[ word ] ) {
        candidates.insert( candidateSeq );
      }
    });

    SequenceList list;
    for( auto &candy : candidates ) {
      list.push_back( *candy );
    }

    return list;
  }

  size_t NumWords() const {
    return mWordToSequences.size();
  }

private:

  void ForEveryWordInSequence( const Sequence &seq, const std::function< void( std::string& ) > &block ) {
    int wordLength = std::min< int >( seq.Length(), mWordSize );
    for( int i = 0; i <= seq.Length() - wordLength; i++ ) {
      std::string word = seq.sequence.substr( i, wordLength );
      ForEveryWordVariation( block, word );
    }
  }

  void ForEveryWordVariation( const std::function< void( std::string& ) > &block, std::string &word, size_t start = 0 )
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
            ForEveryWordVariation( block, word, pos + 1 );
          }
        }
      }
    } else {
      block( word );
    }
  }

  int mWordSize;

  SequenceList mSequences;
  std::unordered_map< std::string, std::deque< const Sequence* >, SequenceHasher > mWordToSequences;
};
