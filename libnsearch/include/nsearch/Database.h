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

  using SequenceNo = uint32_t; // SequenceNo

  std::vector< uint32_t > mIndexByKmer;
  std::vector< uint32_t > mNumEntriesByKmer;
  std::vector< SequenceNo > mSequenceNoByKmer;

  OnProgressCallback mProgressCallback;
};
