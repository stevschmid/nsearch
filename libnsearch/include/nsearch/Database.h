#pragma once

#include <unordered_map>
#include <unordered_set>

#include "Sequence.h"
#include "Utils.h"
#include "Aligner.h"
#include "Kmer.h"

#include "Alignment/Seed.h"
#include "Alignment/Ranges.h"
#include "Alignment/HitTracker.h"
#include "Alignment/OptimalChainFinder.h"
#include "Alignment/ExtendAlign.h"

class Database {
  using SequenceRef = std::shared_ptr< Sequence >;
  using SequenceInfo = std::pair< size_t, SequenceRef >;
  using SequenceMappingDatabase = std::unordered_map< Sequence, std::deque< SequenceInfo > >;

public:
  Database( int indexingWordSize, int alignmentSeedSize )
    : mIndexingWordSize( indexingWordSize ), mAlignmentSeedSize( alignmentSeedSize )
  {
    AlignmentParams ap;
    ap.xDrop = 32;
    mDP = ExtendAlign( ap );

    // Sequences have to fit in the calculated hash
    assert( indexingWordSize * 2 <= sizeof( size_t ) * 8 );
    assert( alignmentSeedSize * 2 <= sizeof( size_t ) * 8 );
  }

  void AddSequence( const Sequence &seq ) {
    // Save
    mSequences.push_back( std::make_shared< Sequence >( seq ) );
    SequenceRef ref = mSequences.back();

    {
      // Kmers for Indexing
      Kmers kmers( *ref, mIndexingWordSize );
      kmers.ForEach( [&]( const Sequence &kmer, size_t pos ) {
        this->mIndexingDB[ kmer ].push_back( SequenceInfo( pos, ref ) );
      });
    }

    {
      // Kmers for Seeding
      Kmers kmers( *ref, mAlignmentSeedSize );
      kmers.ForEach( [&]( const Sequence &kmer, size_t pos ) {
        this->mSeedDB[ ref ][ kmer ].push_back( SequenceInfo( pos, ref ) );
      });
    }
  }

  SequenceList Query( const Sequence &query, int maxHits = 10 ) {
    std::unordered_set< SequenceRef > rawCandidates;

    // Go through each kmer, find candidates
    Kmers kmers( query, mIndexingWordSize );
    kmers.ForEach( [&]( const Sequence &kmer, size_t pos ) {
      for( auto &seqInfo : mIndexingDB[ kmer ] ) {
        SequenceRef seq = seqInfo.second;
        rawCandidates.insert( seq );
      }
    });

    // For each candidate, find alignment
    std::unordered_map< SequenceRef, HitTracker > candidates;
    Kmers kmers2( query, mAlignmentSeedSize );
    kmers2.ForEach( [&]( const Sequence &kmer, size_t pos ) {
      for( auto &candidate : rawCandidates ) {
        for( auto &seqInfo : mSeedDB[ candidate ][ kmer ] ) {
          size_t cpos = seqInfo.first;
          SequenceRef cseq = seqInfo.second;
          candidates[ cseq ].AddHit( pos, cpos, mAlignmentSeedSize );
        }
      }
    });

    for( auto &c : candidates ) {
      const Sequence &candidate = *c.first;
      const HitTracker &hitTracker = c.second;

      for( auto &seed : hitTracker.Seeds() ) {
        size_t leftQuery, leftCandidate;
        int left = mDP.Extend( query, candidate,
            ExtendAlign::ExtendDirection::backwards,
            seed.s1, seed.s2,
            &leftQuery, &leftCandidate);

        size_t rightQuery, rightCandidate;
        int right = mDP.Extend( query, candidate,
            ExtendAlign::ExtendDirection::forwards,
            seed.s1 + seed.length, seed.s2 + seed.length,
            &rightQuery, &rightCandidate );

        Seed ext = seed;
        ext.s1 = leftQuery;
        ext.s2 = leftCandidate;
        ext.length = rightQuery - ext.s1;
        if( ext.length != seed.length ) {
          std::cout << query.Subsequence( ext.s1, ext.length ) << std::endl;
          std::cout << candidate.Subsequence( ext.s2, ext.length ) << std::endl;
          std::cout << "Length: " << ext.length << std::endl;
          std::cout << "=========" << std::endl;
        }

      }
    }
    return SequenceList();

    // Extend each seed -> HSP

    // Try to join HSPs together

    // Sort candidates based on the score
    /* for( auto &c : candidates ) { */
    /*   SequenceRef seq = c.first; */
    /*   const HitTracker &hitTracker = c.second; */

    /*   // Compute optimal chain */
    /*   OptimalChainFinder ocf( hitTracker.Seeds() ); */
    /*   highscore.insert( std::pair< OptimalChainFinder, SequenceRef >( ocf, seq ) ); */
    /* } */

    /* Alignment aln; */
    /* SequenceList list; */

    /* for( auto it = highscore.rbegin(); it != highscore.rend(); ++it ) { */
    /*   const OptimalChainFinder &ocf = it->first; */
    /*   const Sequence &reference = *it->second; */

    /*   int score = mDP.AlignAlongChain( query, reference, ocf.OptimalChain(), &aln ); */
    /*   list.push_back( *(*it).second ); */

    /*   std::cout << "Query " << query.identifier << std::endl; */
    /*   std::cout << "Reference " << reference.identifier << std::endl; */
    /*   std::cout << aln << std::endl; */
    /*   std::cout << " Chain Score " << ocf.Score() */
    /*     << std::endl << " Align Score: " << score */
    /*     << std::endl << " Ref Length " << reference.Length() << std::endl; */
    /*   /1* mDP.DebugPrint( true ); *1/ */

    /*   if( list.size() >= maxHits ) */
    /*     break; */
    /* } */

    /* return list; */
  }

private:
  /* XDropExtension mDP; */
  ExtendAlign mDP;
  size_t mIndexingWordSize, mAlignmentSeedSize;
  std::deque< SequenceRef > mSequences;

  SequenceMappingDatabase mIndexingDB;
  std::unordered_map< SequenceRef, SequenceMappingDatabase > mSeedDB;
};
