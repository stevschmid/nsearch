#include "nsearch/Database.h"

Database::Database( const SequenceList &sequences, size_t kmerLength, const OnProgressCallback &progressCallback )
  : mSequences( sequences ), mKmerLength( kmerLength ), mProgressCallback( progressCallback )
{
  mMaxUniqueKmers = 1 << ( 2 * mKmerLength ); // 2 bits per nt

  size_t totalEntries = 0;
  size_t totalUniqueEntries = 0;

  std::vector< uint32_t > uniqueCount( mMaxUniqueKmers );
  std::vector< SequenceNo > uniqueIndex( mMaxUniqueKmers, -1 );

  for( SequenceNo seqNo = 0; seqNo < mSequences.size(); seqNo++ ) {
    const Sequence &seq = mSequences[ seqNo ];

    Kmers kmers( seq, mKmerLength );
    kmers.ForEach( [&]( Kmer kmer, size_t pos ) {
      totalEntries++;

      // Count unique words
      if( uniqueIndex[ kmer ] == seqNo )
        return;

      uniqueIndex[ kmer ] = seqNo;
      uniqueCount[ kmer ]++;
      totalUniqueEntries++;
    });

    // Progress
    if( seqNo % 500 == 0 || seqNo + 1 == mSequences.size() ) {
      mProgressCallback( ProgressType::StatsCollection, seqNo + 1, mSequences.size() );
    }
  } // Calculate indices

  mIndexByKmer.reserve( mMaxUniqueKmers );
  for( size_t i = 0; i < mMaxUniqueKmers; i++ ) {
    mIndexByKmer[ i ] = i > 0 ? mIndexByKmer[ i - 1 ] + uniqueCount[ i - 1 ] : 0;
  }

  // Populate DB
  mSequenceNoByKmer.resize( totalUniqueEntries );

  // Reset to 0
  mNumEntriesByKmer = std::vector< uint32_t >( mMaxUniqueKmers );

  uniqueIndex = std::vector< SequenceNo>( mMaxUniqueKmers, -1 );

  for( SequenceNo seqNo = 0; seqNo < mSequences.size(); seqNo++ ) {
    const Sequence &seq = mSequences[ seqNo ];

    Kmers kmers( seq, mKmerLength );
    kmers.ForEach( [&]( Kmer kmer, size_t pos ) {
      if( uniqueIndex[ kmer ] == seqNo )
        return;

      uniqueIndex[ kmer ] = seqNo;

      mSequenceNoByKmer[ mIndexByKmer[ kmer ] + mNumEntriesByKmer[ kmer ] ] = seqNo;
      mNumEntriesByKmer[ kmer ]++;
    });

    // Progress
    if( seqNo % 500 == 0 || seqNo + 1 == mSequences.size() ) {
      mProgressCallback( ProgressType::Indexing, seqNo + 1, mSequences.size() );
    }
  }
}

size_t Database::Size() const {
  return mSequences.size();
}
