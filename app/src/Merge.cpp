#include "Merge.h"

#include <nsearch/Sequence.h>
#include <nsearch/FASTQ/Writer.h>
#include <nsearch/FASTA/Reader.h>
#include <nsearch/PairedEnd/Merger.h>
#include <nsearch/PairedEnd/Reader.h>

#include "Common.h"
#include "Stats.h"
#include "WorkerQueue.h"

template<>
class QueueItemInfo< SequenceList > {
public:
  static size_t Count( const SequenceList &list ) { return list.size(); }
};

class MergedReadWriterWorker {
public:
  MergedReadWriterWorker( const std::string &path )
    : mWriter( path )
  {
  }

  void Process( const SequenceList &queueItem ) {
    for( auto seq : queueItem ) {
      mWriter << seq;
    }
  }

private:
  FASTQ::Writer mWriter;
};
using MergedReadWriter = WorkerQueue< SequenceList, MergedReadWriterWorker, const std::string& >;

using PairedReads = std::pair< SequenceList, SequenceList >;

template<>
class QueueItemInfo< PairedReads > {
public:
  static size_t Count( const PairedReads &list ) { return list.first.size(); }
};

class ReadMergerWorker {
public:
  ReadMergerWorker( MergedReadWriter& writer )
    : mWriter( writer )
  {
  }

  void Process( const PairedReads &queueItem ) {
    const SequenceList &fwd = queueItem.first;
    const SequenceList &rev = queueItem.second;

    const PairedEnd::Merger &merger = mMerger;

    Sequence mergedRead;
    SequenceList mergedReads;

    auto fit = fwd.begin();
    auto rit = rev.begin();
    while( fit != fwd.end() && rit != rev.end() ) {
      if( merger.Merge( mergedRead, *fit, *rit ) ) {
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
  MergedReadWriter &mWriter;
  PairedEnd::Merger mMerger;
};
using ReadMerger = WorkerQueue< PairedReads, ReadMergerWorker, MergedReadWriter& >;

bool Merge( const std::string &fwdPath, const std::string &revPath, const std::string &mergedPath ) {
  const int numReadsPerWorkItem = 512;

  PairedEnd::Reader reader( fwdPath, revPath );

  MergedReadWriter writer( 1, mergedPath );
  ReadMerger merger( -1, writer );

  SequenceList fwdReads, revReads;

  enum ProgressType { ReadFile, MergeReads, WriteReads };

  ProgressOutput progress;
  progress.Add( ProgressType::ReadFile, "Reading files", UnitType::BYTES );
  progress.Add( ProgressType::MergeReads, "Merging reads" );
  progress.Add( ProgressType::WriteReads, "Writing merged reads" );

  merger.OnProcessed( [&]( size_t numProcessed, size_t numEnqueued ) {
    progress.Set( ProgressType::MergeReads, numProcessed, numEnqueued );
  });

  writer.OnProcessed( [&]( size_t numProcessed, size_t numEnqueued ) {
    progress.Set( ProgressType::WriteReads, numProcessed, numEnqueued );
  });

  progress.Activate( ProgressType::ReadFile );
  while( !reader.EndOfFile() ) {
    reader.Read( fwdReads, revReads, numReadsPerWorkItem );
    auto pair = std::pair< SequenceList, SequenceList >( std::move( fwdReads ), std::move( revReads ) );
    merger.Enqueue( pair );
    progress.Set( ProgressType::ReadFile, reader.NumBytesRead(), reader.NumBytesTotal() );
  }

  progress.Activate( ProgressType::MergeReads );
  merger.WaitTillDone();

  progress.Activate( ProgressType::WriteReads );
  writer.WaitTillDone();

  return true;
}
