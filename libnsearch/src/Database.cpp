#include "nsearch/Database.h"

Database::Database( const SequenceList &sequences, size_t kmerLength, const OnProgressCallback &progressCallback )
  : mSequences( sequences ), mKmerLength( kmerLength ), mProgressCallback( progressCallback )
{
  mMaxUniqueKmers = 1 << ( 2 * mKmerLength ); // 2 bits per nt

  size_t totalEntries = 0;
  size_t totalUniqueEntries = 0;

  std::vector< uint32_t > uniqueCount( mMaxUniqueKmers );
  std::vector< SequenceId > uniqueIndex( mMaxUniqueKmers, -1 );

  for( SequenceId seqId = 0; seqId < mSequences.size(); seqId++ ) {
    const Sequence &seq = mSequences[ seqId ];

    Kmers kmers( seq, mKmerLength );
    kmers.ForEach( [&]( Kmer kmer, size_t pos ) {
      totalEntries++;

      // Count unique words
      if( uniqueIndex[ kmer ] == seqId )
        return;

      uniqueIndex[ kmer ] = seqId;
      uniqueCount[ kmer ]++;
      totalUniqueEntries++;
    });

    // Progress
    if( seqId % 500 == 0 || seqId + 1 == mSequences.size() ) {
      mProgressCallback( ProgressType::StatsCollection, seqId + 1, mSequences.size() );
    }
  } // Calculate indices

  mSequenceIdsOffsetByKmer.reserve( mMaxUniqueKmers );
  for( size_t i = 0; i < mMaxUniqueKmers; i++ ) {
    mSequenceIdsOffsetByKmer[ i ] = i > 0 ? mSequenceIdsOffsetByKmer[ i - 1 ] + uniqueCount[ i - 1 ] : 0;
  }

  // Populate DB
  mSequenceIds.reserve( totalUniqueEntries );

  // Reset to 0
  mSequenceIdsCountByKmer = std::vector< uint32_t >( mMaxUniqueKmers );

  uniqueIndex = std::vector< SequenceId>( mMaxUniqueKmers, -1 );

  for( SequenceId seqId = 0; seqId < mSequences.size(); seqId++ ) {
    const Sequence &seq = mSequences[ seqId ];

    Kmers kmers( seq, mKmerLength );
    kmers.ForEach( [&]( Kmer kmer, size_t pos ) {
      if( uniqueIndex[ kmer ] == seqId )
        return;

      uniqueIndex[ kmer ] = seqId;

      mSequenceIds[ mSequenceIdsOffsetByKmer[ kmer ] + mSequenceIdsCountByKmer[ kmer ] ] = seqId;
      mSequenceIdsCountByKmer[ kmer ]++;
    });

    // Progress
    if( seqId % 500 == 0 || seqId + 1 == mSequences.size() ) {
      mProgressCallback( ProgressType::Indexing, seqId + 1, mSequences.size() );
    }
  }
}

size_t Database::Size() const {
  return mSequences.size();
}
