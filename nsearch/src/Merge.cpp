#include "Merge.h"

#include <nsearch/PairedEnd/Merger.h>
#include <nsearch/PairedEnd/Reader.h>
#include <nsearch/Sequence.h>
#include <nsearch/Alphabet/DNA.h>

#include <memory>

#include "Common.h"
#include "FileFormat.h"
#include "Stats.h"
#include "WorkerQueue.h"

template < typename A >
class QueueItemInfo< SequenceList< A > > {
public:
  static size_t Count( const SequenceList< A >& list ) {
    return list.size();
  }
};

template < typename A >
class MergedReadWriterWorker {
public:
  MergedReadWriterWorker( const std::string& path )
      : mWriter( std::move(
          DetectFileFormatAndOpenWriter< A >( path, FileFormat::FASTQ ) ) ) {}

  void Process( const SequenceList< A >& queueItem ) {
    for( auto seq : queueItem ) {
      ( *mWriter ) << seq;
    }
  }

private:
  std::unique_ptr< SequenceWriter< A > > mWriter;
};

template < typename A >
using MergedReadWriter = WorkerQueue< MergedReadWriterWorker< A >,
                                      SequenceList< A >, const std::string& >;

template < typename A >
using PairedReads = std::pair< SequenceList< A >, SequenceList< A > >;

template < typename A >
class QueueItemInfo< PairedReads< A > > {
public:
  static size_t Count( const PairedReads< A >& list ) {
    return list.first.size();
  }
};

template < typename A >
class ReadMergerWorker {
public:
  ReadMergerWorker( MergedReadWriter< A >* writer ) : mWriter( *writer ) {}

  void Process( const PairedReads< A >& queueItem ) {
    const SequenceList< A >& fwd = queueItem.first;
    const SequenceList< A >& rev = queueItem.second;

    const PairedEnd::Merger< A >& merger = mMerger;

    Sequence< A >     mergedRead;
    SequenceList< A > mergedReads;

    auto fit = fwd.begin();
    auto rit = rev.begin();
    while( fit != fwd.end() && rit != rev.end() ) {
      if( merger.Merge( *fit, *rit, &mergedRead ) ) {
        gStats.numMerged++;
        gStats.mergedReadsTotalLength += mergedRead.Length();
        mergedReads.push_back( std::move( mergedRead ) );
      }

      ++fit;
      ++rit;
    }

    if( !mergedReads.empty() )
      mWriter.Enqueue( mergedReads );

    gStats.numProcessed += fwd.size();
  }

private:
  MergedReadWriter< A >& mWriter;
  PairedEnd::Merger< A > mMerger;
};

template < typename A >
using ReadMerger = WorkerQueue< ReadMergerWorker< A >, PairedReads< A >,
                                MergedReadWriter< A >* >;

bool DoMerge( const std::string& fwdPath, const std::string& revPath,
              const std::string& mergedPath ) {
  const int numReadsPerWorkItem = 512;

  PairedEnd::Reader< DNA > reader( fwdPath, revPath );
  MergedReadWriter< DNA >  writer( 1, mergedPath );
  ReadMerger< DNA >        merger( -1, &writer );

  SequenceList< DNA > fwdReads, revReads;

  enum ProgressType { ReadFile, MergeReads, WriteReads };

  ProgressOutput progress;
  progress.Add( ProgressType::ReadFile, "Read files", UnitType::BYTES );
  progress.Add( ProgressType::MergeReads, "Merge reads" );
  progress.Add( ProgressType::WriteReads, "Write merged reads" );

  merger.OnProcessed(
    [&]( const size_t numProcessed, const size_t numEnqueued ) {
      progress.Set( ProgressType::MergeReads, numProcessed, numEnqueued );
    } );

  writer.OnProcessed(
    [&]( const size_t numProcessed, const size_t numEnqueued ) {
      progress.Set( ProgressType::WriteReads, numProcessed, numEnqueued );
    } );

  progress.Activate( ProgressType::ReadFile );
  while( !reader.EndOfFile() ) {
    reader.Read( numReadsPerWorkItem, &fwdReads, &revReads );
    auto pair = std::pair< SequenceList< DNA >, SequenceList< DNA > >(
      std::move( fwdReads ), std::move( revReads ) );
    merger.Enqueue( pair );
    progress.Set( ProgressType::ReadFile, reader.NumBytesRead(),
                  reader.NumBytesTotal() );
  }

  progress.Activate( ProgressType::MergeReads );
  merger.WaitTillDone();

  progress.Activate( ProgressType::WriteReads );
  writer.WaitTillDone();

  return true;
}
