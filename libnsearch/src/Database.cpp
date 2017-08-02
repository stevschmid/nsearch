#include "nsearch/Database.h"

Database::Database( const size_t kmerLength )
    : mKmerLength( kmerLength ),
      mProgressCallback( []( ProgressType, const size_t, const size_t ) {} ),
      mMaxUniqueKmers( 1 << ( 2 * mKmerLength ) ) // 2 bits per nt
{}

void Database::SetProgressCallback(
  const OnProgressCallback& progressCallback ) {
  mProgressCallback = progressCallback;
}

void Database::Initialize( const SequenceList& sequences ) {
  mSequences = sequences;

  size_t totalEntries       = 0;
  size_t totalUniqueEntries = 0;

  /* std::vector< uint32_t > count( mMaxUniqueKmers ); */
  std::vector< size_t >     uniqueCount( mMaxUniqueKmers );
  std::vector< SequenceId > uniqueIndex( mMaxUniqueKmers, -1 );

  for( SequenceId seqId = 0; seqId < mSequences.size(); seqId++ ) {
    const Sequence& seq = mSequences[ seqId ];

    Kmers kmers( seq, mKmerLength );
    kmers.ForEach( [&]( const Kmer kmer, const size_t pos ) {
      totalEntries++;

      // Count unique words
      if( uniqueIndex[ kmer ] == seqId )
        return;

      uniqueIndex[ kmer ] = seqId;
      uniqueCount[ kmer ]++;
      totalUniqueEntries++;
    } );

    // Progress
    if( seqId % 512 == 0 || seqId + 1 == mSequences.size() ) {
      mProgressCallback( ProgressType::StatsCollection, seqId + 1,
                         mSequences.size() );
    }
  } // Calculate indices

  mSequenceIdsOffsetByKmer.reserve( mMaxUniqueKmers );
  for( size_t i = 0; i < mMaxUniqueKmers; i++ ) {
    mSequenceIdsOffsetByKmer[ i ] =
      i > 0 ? mSequenceIdsOffsetByKmer[ i - 1 ] + uniqueCount[ i - 1 ] : 0;
  }

  // Populate DB
  mSequenceIds.reserve( totalUniqueEntries );
  mKmers.reserve( totalEntries );

  // Reset to 0
  mSequenceIdsCountByKmer = std::vector< size_t >( mMaxUniqueKmers );
  mKmerCountBySequenceId  = std::vector< size_t >( mSequences.size() );
  mKmerOffsetBySequenceId = std::vector< size_t >( mSequences.size() );

  uniqueIndex = std::vector< SequenceId >( mMaxUniqueKmers, -1 );

  auto   kmersData = mKmers.data();
  size_t kmerCount = 0;

  for( SequenceId seqId = 0; seqId < mSequences.size(); seqId++ ) {
    const Sequence& seq = mSequences[ seqId ];

    mKmerOffsetBySequenceId[ seqId ] = kmerCount;

    Kmers kmers( seq, mKmerLength );
    kmers.ForEach( [&]( const Kmer kmer, const size_t pos ) {
      // Encode position in kmersData implicitly
      // by saving _every_ kmer
      kmersData[ kmerCount++ ] = kmer;

      if( uniqueIndex[ kmer ] == seqId )
        return;

      uniqueIndex[ kmer ] = seqId;

      mSequenceIds[ mSequenceIdsOffsetByKmer[ kmer ] +
                    mSequenceIdsCountByKmer[ kmer ] ] = seqId;
      mSequenceIdsCountByKmer[ kmer ]++;
    } );

    mKmerCountBySequenceId[ seqId ] =
      kmerCount - mKmerOffsetBySequenceId[ seqId ];

    // Progress
    if( seqId % 512 == 0 || seqId + 1 == mSequences.size() ) {
      mProgressCallback( ProgressType::Indexing, seqId + 1, mSequences.size() );
    }
  }
}

const Sequence& Database::GetSequenceById( const SequenceId& seqId ) const {
  assert( seqId < NumSequences() );
  return mSequences[ seqId ];
}

size_t Database::NumSequences() const {
  return mSequences.size();
}

size_t Database::MaxUniqueKmers() const {
  return mMaxUniqueKmers;
}

size_t Database::KmerLength() const {
  return mKmerLength;
}

bool Database::GetKmersForSequenceId( const SequenceId& seqId,
                                      const Kmer**      kmers,
                                      size_t*           numKmers ) const {
  if( seqId >= NumSequences() )
    return false;

  const auto& offset = mKmerOffsetBySequenceId[ seqId ];
  const auto& count  = mKmerCountBySequenceId[ seqId ];

  *kmers    = &mKmers[ offset ];
  *numKmers = count;
  return count > 0;
}

bool Database::GetSequenceIdsIncludingKmer( const Kmer&        kmer,
                                            const SequenceId** seqIds,
                                            size_t* numSeqIds ) const {
  if( kmer >= MaxUniqueKmers() )
    return false;

  const auto& offset = mSequenceIdsOffsetByKmer[ kmer ];
  const auto& count  = mSequenceIdsCountByKmer[ kmer ];

  *seqIds    = &mSequenceIds[ offset ];
  *numSeqIds = count;
  return count > 0;
}
