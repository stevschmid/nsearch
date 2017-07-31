#pragma once

#include <deque>
#include <vector>

#include "Sequence.h"
#include "Utils.h"

#include "Database/HSP.h"
#include "Database/Highscore.h"
#include "Database/Kmers.h"

class Database {
  friend class GlobalSearch;

public:
  enum ProgressType { StatsCollection, Indexing };
  using OnProgressCallback =
    std::function< void( ProgressType, const size_t, const size_t ) >;

  Database( const SequenceList& sequences, const size_t kmerLength,
            const OnProgressCallback& progressCallback =
              []( ProgressType, const size_t, const size_t ) {} );
  size_t Size() const;

private:
  size_t mKmerLength;

  SequenceList mSequences;
  size_t       mMaxUniqueKmers;

  using SequenceId = uint32_t; // SequenceId

  std::vector< size_t >     mSequenceIdsOffsetByKmer;
  std::vector< size_t >     mSequenceIdsCountByKmer;
  std::vector< SequenceId > mSequenceIds;

  std::vector< size_t > mKmerOffsetBySequenceId;
  std::vector< size_t > mKmerCountBySequenceId;
  std::vector< Kmer >   mKmers;

  OnProgressCallback mProgressCallback;
};
