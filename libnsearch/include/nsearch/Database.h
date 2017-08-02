#pragma once

#include <deque>
#include <vector>

#include "Sequence.h"
#include "Utils.h"

#include "Database/HSP.h"
#include "Database/Highscore.h"
#include "Database/Kmers.h"

using SequenceId = uint32_t; // SequenceId

template < typename Alphabet >
class Database {
public:
  enum ProgressType { StatsCollection, Indexing };
  using OnProgressCallback =
    std::function< void( ProgressType, const size_t, const size_t ) >;

  Database( const size_t kmerLength );

  void SetProgressCallback( const OnProgressCallback& progressCallback );
  void Initialize( const SequenceList< Alphabet >& sequences ) const;

  size_t NumSequences() const;
  size_t MaxUniqueKmers() const;
  size_t KmerLength() const;

  const Sequence< Alphabet >& GetSequenceById( const SequenceId& seqId ) const;

  bool GetKmersForSequenceId( const SequenceId& seqId, const Kmer** kmers,
                              size_t* numKmers ) const;
  bool GetSequenceIdsIncludingKmer( const Kmer& kmer, const SequenceId** seqIds,
                                    size_t* numSeqIds ) const;

private:
  size_t mKmerLength;

  SequenceList< Alphabet > mSequences;
  size_t                   mMaxUniqueKmers;

  std::vector< size_t >     mSequenceIdsOffsetByKmer;
  std::vector< size_t >     mSequenceIdsCountByKmer;
  std::vector< SequenceId > mSequenceIds;

  std::vector< size_t > mKmerOffsetBySequenceId;
  std::vector< size_t > mKmerCountBySequenceId;
  std::vector< Kmer >   mKmers;

  OnProgressCallback mProgressCallback;
};

#include "../../src/Database.cpp"

