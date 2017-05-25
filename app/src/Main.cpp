#include <docopt.h>
#include <iostream>

#include <nsearch/Sequence.h>
#include <nsearch/FASTQ/Writer.h>
#include <nsearch/PairedEnd/Merger.h>
#include <nsearch/PairedEnd/Reader.h>

#include "Stats.h"
#include "ThreadPool.h"

Stats gStats;

static const char USAGE[] =
R"(
  Process and search sequences.

  Usage:
    nsearch merge <forward.fastq> <reverse.fastq> <merged.fastq>
)";

void PrintProgressLine( size_t numProcessedReads, size_t numTotalReads )
{
  static int PROGRESS_BAR_WIDTH = 50;

  double done = double( numProcessedReads ) / double( numTotalReads );
  int pos = done * PROGRESS_BAR_WIDTH;

  std::cout << "\r[";
  std::cout << std::string( pos, '=' );
  if( done < 1.0 )
    std::cout << '>';
  std::cout << std::string( std::max( PROGRESS_BAR_WIDTH - pos - 1, 0 ), ' ' );
  std::cout << "] ";
  printf( " %.1f%%", done * 100.0 );
  std::cout << std::flush;
}

void PrintSummaryLine( double value, const char *line, double total = 0.0 )
{
  printf( "%10.1f  %s", value, line );
  if( total > 0.0 ) {
    printf( " (%.1f%%)", value / total * 100.0 );
  }
  printf( "\n" );
}

class ThreadSafeWriter {
private:
  std::queue< SequenceList > mQueue;
  std::mutex mQueueMutex;
  FASTQ::Writer mWriter;
  std::atomic< bool > mStop;

  std::atomic< bool > mWorking;
  std::thread mThread;

  void Loop() {
    SequenceList list;

    while( !mStop )
    {
      { // acquire lock
        std::unique_lock< std::mutex > lock( mQueueMutex );
        if( !mQueue.empty() ) {
          list = std::move( mQueue.front() );
          mQueue.pop();
        }
        mWorking = true;
      } // release lock

      for( Sequence &seq : list ) {
        mWriter << seq;
      }
      list.clear();

      mWorking = false;
    }

    std::this_thread::sleep_for( std::chrono::milliseconds( 50 ) );
  }

public:
  ThreadSafeWriter( const std::string &path )
    : mWriter( path ), mStop( false ), mWorking( false )
  {
    mThread = std::thread( [this] { this->Loop(); } );
  }

  ~ThreadSafeWriter() {
    mStop = true;
    if( mThread.joinable() )
      mThread.join();
  }

  void Enqueue( SequenceList &list ) {
    std::lock_guard< std::mutex > lock( mQueueMutex );
    mQueue.push( std::move( list ) );
  }

  bool Done() const {
    return !mWorking && mQueue.empty();
  }
};

bool Merge( const std::string &fwdPath, const std::string &revPath, const std::string &mergedPath ) {
  ThreadSafeWriter writer( mergedPath );

  PairedEnd::Merger mergerObj;
  const PairedEnd::Merger &merger = mergerObj; // const so we ensure thread-safety

  PairedEnd::Reader reader( fwdPath, revPath );

  ThreadPool pool;
  SequenceList fwdReads, revReads;

  while( !reader.EndOfFile() ) {
    reader.Read( fwdReads, revReads, 512 );

    auto mergeAndWrite = [ &writer, &merger ] ( SequenceList &fwd, SequenceList &rev ) {
      SequenceList mergedReads;
      Sequence mergedRead;

      auto fit = fwd.begin();
      auto rit = rev.begin();
      while( fit != fwd.end() && rit != rev.end() ) {
        if( merger.Merge( mergedRead, *fit, *rit ) ) {
          gStats.numMerged++;
          gStats.mergedReadsTotalLength += mergedRead.Length();

          mergedReads.push_back( std::move( mergedRead ) );
        }

        gStats.numProcessed++;

        ++fit;
        ++rit;
      }

      if( !mergedReads.empty() )
        writer.Enqueue( mergedReads );
    };

    auto task = std::bind( mergeAndWrite, std::move( fwdReads ), std::move( revReads ) );
    pool.Enqueue( task );
    PrintProgressLine( reader.NumBytesRead(), reader.NumBytesTotal() );
  }

  while( !writer.Done() && !pool.Done() )
    std::this_thread::sleep_for( std::chrono::milliseconds( 50 ) );

  return true;
}

int main( int argc, const char **argv ) {
  std::map<std::string, docopt::value> args
    = docopt::docopt(USAGE,
        { argv + 1, argv + argc },
        true, // help
        "nsearch");

  if( args["merge"] ) {
    gStats.StartTimer();

    Merge( args[ "<forward.fastq>" ].asString(),
        args[ "<reverse.fastq>" ].asString(),
        args[ "<merged.fastq>" ].asString() );


    gStats.StopTimer();

    std::cout << std::endl;
    std::cout << "Summary:" << std::endl;
    PrintSummaryLine( gStats.ElapsedMillis() / 1000.0, "Seconds" );
    PrintSummaryLine( gStats.numProcessed / gStats.ElapsedMillis(), "Processed/ms" );
    PrintSummaryLine( gStats.numProcessed, "Pairs" );
    PrintSummaryLine( gStats.numMerged, "Merged", gStats.numProcessed );
    PrintSummaryLine( gStats.MeanMergedLength(), "Mean merged length" );
  }

  return 0;
}
