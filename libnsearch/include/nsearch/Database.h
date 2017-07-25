#pragma once

#include <vector>
#include <deque>

#include "Sequence.h"
#include "Utils.h"

#include "Database/Kmers.h"
#include "Database/HSP.h"
#include "Database/Highscore.h"

class Database {
  friend class GlobalSearch;

public:
  enum ProgressType { StatsCollection, Indexing };
  using OnProgressCallback = std::function< void ( ProgressType, size_t, size_t )  >;

  Database( const SequenceList &sequences, size_t wordSize,
      const OnProgressCallback &progressCallback = []( ProgressType, size_t, size_t ) { } );
  size_t Size() const;

private:
  size_t mWordSize;

  SequenceList mSequences;
  size_t mMaxUniqueWords;

  using WordEntry = struct WordEntry_s {
    uint32_t sequence;
    uint32_t pos;
    WordEntry_s *nextEntry;

    WordEntry_s()
    : sequence( -1 )
    {
    }

    WordEntry_s( uint32_t s, uint32_t p, WordEntry_s *ne = NULL )
      : sequence( s ), pos( p ), nextEntry( ne )
    {

    }
  };

  std::vector< uint32_t > mIndexByWord;
  std::vector< uint32_t > mNumEntriesByWord;
  std::vector< WordEntry > mFirstEntries; // first word (kmer) hit for each candidate
  std::vector< WordEntry > mFurtherEntries;

  OnProgressCallback mProgressCallback;
};
