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

  Database( const SequenceList &sequences, size_t kmerLength,
      const OnProgressCallback &progressCallback = []( ProgressType, size_t, size_t ) { } );
  size_t Size() const;

private:
  size_t mKmerLength;

  SequenceList mSequences;
  size_t mMaxUniqueKmers;

  using SequenceId = uint32_t; // SequenceId

  using KmerInfo = struct KmerInfo_s {
    size_t pos;
    Kmer kmer;
    SequenceId seqId;

    KmerInfo_s *prev, *next;
  };

  std::vector< uint32_t > mSequenceIdsOffsetByKmer;
  std::vector< uint32_t > mSequenceIdsCountByKmer;
  std::vector< SequenceId > mSequenceIds;

  std::vector< KmerInfo > mKmerInfos;

  OnProgressCallback mProgressCallback;
};
