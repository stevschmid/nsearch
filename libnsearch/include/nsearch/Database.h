#pragma once

#include <unordered_map>
#include <unordered_set>
#include <set>

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
  Database( size_t wordSize )
    : mWordSize( wordSize )
  {
    // Sequences have to fit in the calculated hash
    assert( mWordSize * 2 <= sizeof( size_t ) * 8 );
  }

  void AddSequence( const Sequence &seq ) {
    // Save
    mSequences.push_back( std::make_shared< Sequence >( seq ) );
    SequenceRef ref = mSequences.back();

    // Kmers for Indexing
    Kmers kmers( *ref, mWordSize );
    kmers.ForEach( [&]( const Sequence &kmer, size_t pos ) {
      this->mWordDB[ kmer ].push_back( SequenceInfo( pos, ref ) );
    });
  }

  SequenceList Query( const Sequence &query, int maxHits = 10 ) {
    const int xDrop = 32;

    std::unordered_map< SequenceRef, HitTracker > hits;

    // Go through each kmer, find candidates
    Kmers kmers( query, mWordSize );
    kmers.ForEach( [&]( const Sequence &kmer, size_t pos ) {
      for( auto &seqInfo : mWordDB[ kmer ] ) {
        size_t candidatePos = seqInfo.first;
        SequenceRef candidateRef = seqInfo.second;

        hits[ candidateRef ].AddHit( pos, candidatePos, mWordSize );
        /* hits[ candidateRef ].push_back( Seed( pos, candidatePos, mWordSize ) ); */
      }
    });

    // Order candidates by number of shared words (hits)
    using Candidate = std::pair< SequenceRef, HitTracker >;
    struct CandidateCompare {
      bool operator() ( const Candidate &left, const Candidate &right ) const {
        return left.second.Score() < right.second.Score();
      }
    };
    using Highscore = std::set< Candidate , CandidateCompare >;
    Highscore highscore;
    for( auto &h : hits ) {
      highscore.insert( std::make_pair( h.first, h.second ) );
    }

    // For each candidate:
    // - Get HSPs,
    // - Check for good HSP (>= similarity threshold)
    // - Join HSP together
    // - Align
    // - Check similarity
    std::cout << "=====" << std::endl;
    for( auto it = highscore.rbegin(); it != highscore.rend(); ++it ) {
      const Candidate &candidate = *it;
      const Sequence &candidateSeq = *candidate.first;

      // some yolo preliminary optimization
      /* if( candidate.second.Score() < 0.5 * query.Length() ) { */
      /*   continue; */
      /* } */

      // Find all HSP
      // Sort by length
      // Try to find best chain
      // Fill space between with banded align

      for( auto &seed : candidate.second.Seeds() ) {
        std::cout << "Seed " << std::endl;
        size_t leftQuery, leftCandidate;
        int leftScore = mXDropExtender.Extend( query, candidateSeq,
            xDrop,
            AlignExtendDirection::backwards,
            seed.s1, seed.s2,
            &leftQuery, &leftCandidate );

        size_t rightQuery, rightCandidate;
        int rightScore = mXDropExtender.Extend( query, candidateSeq,
            xDrop,
            AlignExtendDirection::forwards,
            seed.s1 + seed.length, seed.s2 + seed.length,
            &rightQuery, &rightCandidate );

        /* ht.AddHit( leftQuery, leftCandidate, rightQuery - leftQuery ); */

        /* int score = seed.length * mDP.AP().matchScore; // initial seed is an exact match */
        /* score += leftScore; // left extended */
        /* score += rightScore; // right extended */

        /* if( ext.length != seed.length ) { */
        /*   std::cout << "Length: " << ext.length << " (ext " */
        /*     << ( ext.length - seed.length ) << ")" << std::endl; */
        /*   std::cout << ext.s1 << " " << ext.s2 << std::endl; */
        /*   /1* std::cout << query.Subsequence( ext.s1, ext.length ) << std::endl; *1/ */
        /*   /1* std::cout << candidate.Subsequence( ext.s2, ext.length ) << std::endl; *1/ */
        /*   std::cout << "=========" << std::endl; */
        /* } */
      }


    }
    return SequenceList();

    /* for( auto &c : candidates ) { */
    /*   SequenceRef cref = c.first; */
    /*   const Sequence &candidate = *cref; */
    /*   const HitTracker &hitTracker = c.second; */

    /*   for( auto &seed : hitTracker.Seeds() ) { */
    /*     size_t leftQuery, leftCandidate; */
    /*     int leftScore = mXDropExtender.Extend( query, candidate, */
    /*         xDrop, */
    /*         AlignExtendDirection::backwards, */
    /*         seed.s1, seed.s2, */
    /*         &leftQuery, &leftCandidate); */

    /*     size_t rightQuery, rightCandidate; */
    /*     int rightScore = mXDropExtender.Extend( query, candidate, */
    /*         xDrop, */
    /*         AlignExtendDirection::forwards, */
    /*         seed.s1 + seed.length, seed.s2 + seed.length, */
    /*         &rightQuery, &rightCandidate ); */

    /*     Seed ext = seed; */
    /*     ext.s1 = leftQuery; */
    /*     ext.s2 = leftCandidate; */
    /*     ext.length = rightQuery - ext.s1; */

    /*     /1* int score = seed.length * mDP.AP().matchScore; // initial seed is an exact match *1/ */
    /*     /1* score += leftScore; // left extended *1/ */
    /*     /1* score += rightScore; // right extended *1/ */

    /*     /1* if( ext.length != seed.length ) { *1/ */
    /*     /1*   std::cout << "Length: " << ext.length << " (ext " *1/ */
    /*     /1*     << ( ext.length - seed.length ) << ")" << std::endl; *1/ */
    /*     /1*   std::cout << ext.s1 << " " << ext.s2 << std::endl; *1/ */
    /*     /1*   /2* std::cout << query.Subsequence( ext.s1, ext.length ) << std::endl; *2/ *1/ */
    /*     /1*   /2* std::cout << candidate.Subsequence( ext.s2, ext.length ) << std::endl; *2/ *1/ */
    /*     /1*   std::cout << "=========" << std::endl; *1/ */
    /*     /1* } *1/ */

    /*   } */
    /* } */

    /* return SequenceList(); */

    // Extend each seed -> HSP

    // Try to join HSPs together

    // Sort candidates based on the score
    /* for( auto &c : candidates ) { */
    /*   SequenceRef seq = c.first; */
    /*   const HitTracker &hitTracker = c.second; */

    /*   // Compute optimal chain */
    /*   OptimalChainFinder ocf( hitTracker.Seeds() ); */
    /*   highscore.insert( std::make_pair( ocf, seq ) ); */
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
  XDropExtendAlign mXDropExtender;
  size_t mWordSize;
  std::deque< SequenceRef > mSequences;

  SequenceMappingDatabase mWordDB;
};
