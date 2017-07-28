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

  /* using KmerInfo = struct KmerInfo_s { */
  /*   uint32_t pos; */
  /*   /1* Kmer kmer; *1/ */
  /*   SequenceId seqId; */

  /*   KmerInfo_s *prev; */
  /*   KmerInfo_s *next; */
  /* }; */

  std::vector< size_t > mSequenceIdsOffsetByKmer;
  std::vector< size_t > mSequenceIdsCountByKmer;
  std::vector< SequenceId > mSequenceIds;

  /* std::vector< uint32_t > mKmerInfosOffsetByKmer; */
  /* std::vector< uint32_t > mKmerInfosCountByKmer; */
  /* std::vector< KmerInfo > mKmerInfos; */

  /* std::vector< uint32_t > mKmerOffsetBySequenceId; */
  /* std::vector< uint32_t > mKmerInfosCountByKmer; */
  /* std::vector< KmerInfo > mKmerInfos; */

  std::vector< size_t > mKmerOffsetBySequenceId;
  std::vector< size_t > mKmerCountBySequenceId;
  std::vector< Kmer > mKmers;

  OnProgressCallback mProgressCallback;
};
