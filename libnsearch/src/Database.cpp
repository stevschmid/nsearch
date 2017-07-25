#include "nsearch/Database.h"

Database::Database( const SequenceList &sequences, size_t wordSize )
  : mSequences( sequences ), mWordSize( wordSize )
{
  mMaxUniqueWords = 1 << ( 2 * mWordSize ); // 2 bits per nt

  size_t totalEntries = 0;
  size_t totalFirstEntries = 0;

  std::vector< uint32_t > uniqueCount( mMaxUniqueWords );
  std::vector< uint32_t > uniqueIndex( mMaxUniqueWords, -1 );
  for( uint32_t idx = 0; idx < mSequences.size(); idx++ ) {
    const Sequence &seq = mSequences[ idx ];

    Kmers kmers( seq, mWordSize );
    kmers.ForEach( [&]( Kmer word, size_t pos ) {
        totalEntries++;

        if( uniqueIndex[ word ] != idx ) {
        uniqueIndex[ word ] = idx;
        uniqueCount[ word ]++;
        totalFirstEntries++;
        }
        });
  } // Calculate indices
  mIndexByWord.reserve( mMaxUniqueWords );
  for( size_t i = 0; i < mMaxUniqueWords; i++ ) {
    mIndexByWord[ i ] = i > 0 ? mIndexByWord[ i - 1 ] + uniqueCount[ i - 1 ] : 0;
  }

  // Populate DB
  mFirstEntries.resize( totalFirstEntries );
  mFurtherEntries.reserve( totalEntries - totalFirstEntries );

  // Reset to 0
  mNumEntriesByWord = std::vector< uint32_t >( mMaxUniqueWords );
  for( uint32_t idx = 0; idx < mSequences.size(); idx++ ) {
    const Sequence &seq = mSequences[ idx ];

    Kmers kmers( seq, mWordSize );
    kmers.ForEach( [&]( Kmer word, size_t pos ) {
        if( mNumEntriesByWord[ word ] == 0 ) {
        // Create new entry
        WordEntry *entry = &mFirstEntries[ mIndexByWord[ word ] ];
        entry->sequence = idx;
        entry->pos = pos;
        entry->nextEntry = NULL;
        mNumEntriesByWord[ word ]++;
        } else {
        // Check if last entry == index
        WordEntry *entry = &mFirstEntries[ mIndexByWord[ word ] + ( mNumEntriesByWord[ word ] - 1 ) ];

        if( entry->sequence == idx ) {
        mFurtherEntries.emplace_back( idx, pos );
        WordEntry *s = entry;
        while( s->nextEntry ) {
        s = s->nextEntry;
        }
        s->nextEntry = &mFurtherEntries.back();
        } else {
        entry++;
        entry->sequence = idx;
        entry->pos = pos;
        entry->nextEntry = NULL;
        mNumEntriesByWord[ word ]++;
        }
        }
    });

  }
}

size_t Database::Size() const {
  return mSequences.size();
}

