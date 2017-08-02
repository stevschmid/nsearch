#include <catch.hpp>

#include <nsearch/Database.h>

#include "Support.h"

TEST_CASE( "Database" ) {
  SequenceList sequences = { "ATGGG", "CATGGCCC", "GAGAGA" };
  Database db( 4 );
  db.Initialize( sequences );

  SECTION( "Sequences" ) {
    REQUIRE( db.NumSequences() == 3 );
    REQUIRE( db.GetSequenceById( 0 ) == Sequence( "ATGGG" ) );
  }

  SECTION( "Kmer support" ) {
    SECTION( "Stats" ) {
      REQUIRE( db.KmerLength() == 4 );
      REQUIRE( db.MaxUniqueKmers() == 256 ); // 4^4
    }

    SECTION( "Search" ) {
      const Database::SequenceId *seqIds;
      size_t numSeqIds;
      bool found;

      found = db.GetSequenceIdsIncludingKmer( Kmerify( "TATA" ), &seqIds,
                                              &numSeqIds );
      REQUIRE( found == false );
      REQUIRE( numSeqIds == 0 );

      found = db.GetSequenceIdsIncludingKmer( Kmerify( "GAGA" ), &seqIds,
                                              &numSeqIds );
      REQUIRE( found == true );
      REQUIRE( numSeqIds == 1 );
      REQUIRE( seqIds[ 0 ] == 2 );

      found = db.GetSequenceIdsIncludingKmer( Kmerify( "ATGG" ), &seqIds,
                                              &numSeqIds );
      REQUIRE( found == true );
      REQUIRE( numSeqIds == 2 );
      REQUIRE( seqIds[ 0 ] == 0 );
      REQUIRE( seqIds[ 1 ] == 1 );
    }
  }

  SECTION( "Kmers for each sequence" ) {
    const Kmer *kmers;
    size_t numKmers;
    bool found;

    found = db.GetKmersForSequenceId( 10, &kmers, &numKmers );
    REQUIRE( found == false );

    found = db.GetKmersForSequenceId( 0, &kmers, &numKmers );
    REQUIRE( found == true );
    REQUIRE( numKmers == 2 );
    REQUIRE( kmers[ 0 ] == Kmerify( "ATGG" ) );
    REQUIRE( kmers[ 1 ] == Kmerify( "TGGG" ) );
  }
}
