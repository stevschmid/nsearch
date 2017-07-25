#pragma once
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <memory>

#include "Sequence.h"
#include "Utils.h"

#include "Database/Kmers.h"
#include "Database/HSP.h"
#include "Database/Highscore.h"

class DatabaseSearcher;

class Database {
  friend class DatabaseSearcher;

public:
  Database( const SequenceList &sequences, size_t wordSize );
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
};
