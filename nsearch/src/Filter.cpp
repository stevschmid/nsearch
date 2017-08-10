#include "Filter.h"

#include <nsearch/Alphabet/DNA.h>
#include <nsearch/Sequence.h>

#include "Common.h"
#include "FileFormat.h"

bool DoFilter( const std::string& inputPath, const std::string& outputPath,
             const float maxExpectedErrors ) {
  auto reader =
    DetectFileFormatAndOpenReader< DNA >( inputPath, FileFormat::FASTQ );

  enum ProgressType { ReadFile, Filter, WriteFile };

  ProgressOutput progress;
  progress.Add( ProgressType::ReadFile, "Read input file", UnitType::BYTES );
  progress.Add( ProgressType::Filter, "Filter reads" );
  progress.Add( ProgressType::WriteFile, "Write final reads" );

  // Read
  SequenceList< DNA > sequences;
  Sequence< DNA >     seq;
  progress.Activate( ProgressType::ReadFile );
  while( !reader->EndOfFile() ) {
    ( *reader ) >> seq;
    sequences.push_back( std::move( seq ) );
    progress.Set( ProgressType::ReadFile, reader->NumBytesRead(),
                  reader->NumBytesTotal() );
  }

  // Filter
  SequenceList< DNA > filtered;
  progress.Activate( ProgressType::Filter );
  size_t count = 0;
  for( auto& s : sequences ) {
    if( s.NumExpectedErrors() <= maxExpectedErrors ) {
      filtered.push_back( std::move( s ) );
    }
    progress.Set( ProgressType::Filter, ++count, sequences.size() );
  }

  // Write
  progress.Activate( ProgressType::WriteFile );
  count = 0;
  auto writer =
    DetectFileFormatAndOpenWriter< DNA >( outputPath, FileFormat::FASTA );
  for( auto& s : filtered ) {
    ( *writer ) << s;
    progress.Set( ProgressType::WriteFile, ++count, filtered.size() );
  }

  return true;
}
