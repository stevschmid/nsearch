#include <docopt.h>
#include <iostream>

#include <nsearch/Sequence.h>
#include <nsearch/FASTQ/Writer.h>
#include <nsearch/PairedEnd/Merger.h>
#include <nsearch/PairedEnd/Reader.h>

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

bool Merge( const std::string &fwd, const std::string &rev, const std::string &merged ) {
  Sequence fwdRead, revRead, mergedRead;

  PairedEnd::Merger merger;
  PairedEnd::Reader reader( fwd, rev );

  FASTQ::Writer writer( merged );

  while( !reader.EndOfFile() ) {
    reader.Read( fwdRead, revRead );
    if( merger.Merge( mergedRead, fwdRead, revRead ) ) {
      writer << mergedRead;
    }

    PrintProgressLine( reader.NumBytesRead(), reader.NumBytesTotal() );
  }

  return true;
}

int main( int argc, const char **argv ) {
  std::map<std::string, docopt::value> args
    = docopt::docopt(USAGE,
        { argv + 1, argv + argc },
        true, // help
        "nsearch");

  if( args["merge"] ) {
    Merge( args[ "<forward.fastq>" ].asString(),
        args[ "<reverse.fastq>" ].asString(),
        args[ "<merged.fastq>" ].asString() );
  }

  return 0;
}
