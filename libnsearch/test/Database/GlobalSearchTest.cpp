#include <catch.hpp>

#include <nsearch/Alphabet/DNA.h>
#include <nsearch/Database.h>
#include <nsearch/Database/GlobalSearch.h>
#include <nsearch/FASTA/Reader.h>

#include <sstream>
#include <string>
#include <vector>

const char DatabaseContents[] = R"(
>RF00807;mir-314;AFFE01007792.1/82767-82854   42026:Drosophila bipectinata
UCGUAACUUGUGUGGCUUCGAAUGUACCUAGUUGAGGAAAAAUCAGUUUG
GAUUUUGUUACCUCUGGUAUUCGAGCCAAUAAGUUCGG
>RF00807;mir-314;AAFS01000446.1/64778-64866   46245:Drosophila pseudoobscura pseudoobscura
UCGUAACUUGUGUGGCUUCGAAUGUACCUAGUUGAGGAAAACUCCGAAAU
GGAUUUUGUUACCUCUGGUAUUCGAGCCAAUAAGUUCGG
>RF00807;mir-314;AANI01017486.1/342740-342830   7244:Drosophila virilis
UCGUAACUUGUGUGGCUUGAAUGUACCUGGUUGAGGAACGAAUUCAACGU
UUGGAUUUUGUUGCCUUUGGUAUUCGAGCCAAUAAGUUCGG
>RF00807;mir-314;AAPU01011627.1/156896-156990   7230:Drosophila mojavensis
UCGUAACUUGUGUGGCUUCGAAUGUACCUCGUCGAGCGAAAAGCGAAUUC
AUUGUUGGAUUUUGUUGCUCUUGGUAUUCGAGCCAAUAAGUUCGG
>RF00752;mir-14;AAWU01029067.1/5309-5242   7176:Culex quinquefasciatus (southern house mosquito)
UGUGGGAGCGAGAUUAAGGCUUGCUGGUUUCACGUUCGAGUAAAGUCAGU
CUUUUUCUCUCUCCUAUU
>RF00752;mir-14;AAGE02012112.1/774-706   7159:Aedes aegypti (yellow fever mosquito)
UGUGGGAGCGAGAUUAAGGCUUGCUGGUCAUUUAUUACACUCGAAGUCAG
UCUUUUUCUCUCUCCUAUU
>RF00752;mir-14;AANI01011011.1/11101-11163   7244:Drosophila virilis
UGUGGGAGCGAGACGGGGACUCACUGUGCUUUUUAUAUAGUCAGUCUUUU
UCUCUCUCCUAUA
>RF00715;mir-383;AAMC01036319.1/13960-13888   8364:Xenopus (Silurana) tropicalis (western clawed frog)
CUCCUCAGAUCAGAAGGUGAUUGUGGCUUUUAGUAGAUAUUAAGCAGCCA
CAGCACUGCCUGGUCAGAAAGAG
>RF00715;mir-383;AAQR03137803.1/1757-1829   30611:Otolemur garnettii (small-eared galago)
CUCCUCAGAUCAGAAGGUGAUUGUGGCUUUGGGUGCAUGGUUAUAAGCCA
CAGCACUGCCUGGUCAGAAAGAG
>RF00715;mir-383;AFEY01400405.1/5982-5910   9305:Sarcophilus harrisii (Tasmanian devil)
CUCCUCAGAUCAGAAGGUGAUUGUGGCUUUGGGCAGACAUGGAACAGCCA
CAUCACUGGCUGGUCAGAAAGAG
>RF01157;sn1185;DQ789405.1/1-66   6238:Caenorhabditis briggsae
AUCGGUGAUGUGAUAUCCAGUUCUGCUACUGAAGCGUUGUGAAGAUUAAC
UUUCCCCGUCUGAGAU
>RF01157;sn1185;AEHI01092347.1/1008-1073   860376:Caenorhabditis angaria
ACUGAUGAUGUUAACUCCAGUUCUGCUACUGAAUGAAUGUGACGAUAUUC
UUUCCCCGACUGAGGU
>RF01157;sn1185;ABLE03029241.1/4849-4913   281687:Caenorhabditis japonica
AUUGAUGAUGUUCAUCCAGUUCUGCUACUGAAUCAGUGUGAAGAUAUUCU
UUCCCCGACUGAGAU
>RF01885;HSR-omega_1;AFPP01029324.1/15839-15914   1041015:Drosophila rhopaloa
ACCACCUAACCAAGCAAUAUGUAUUUCUUUCUCUAAACUUUAUAGUUGGG
CGUUGAAAGUUGAUAUCGAUCCGUGA
>RF01885;HSR-omega_1;AFFH01007186.1/256640-256716   30033:Drosophila kikkawai
AUCACUUAACCAGCAAUAUGUAUUUCUUUCUCUAAACUUUAUAGUUGGGC
GUUGAAAGUUGAUACGCGAACGUGAAA
>RF01885;HSR-omega_1;AAPU01011178.1/205842-205767   7230:Drosophila mojavensis
ACACGUUAACCAAGCAUUAUGUAUUUCUUUCUCUAAACUUUAUAGUUGGG
CGUUGAAAGUUGAUACGCGAUCGAAC
)";

TEST_CASE( "Global Search" ) {
  Database< DNA > db( 8 );

  SearchParams< DNA > sp;
  sp.maxAccepts = 1;
  sp.maxRejects = 8;
  sp.minIdentity = 0.75f;

  std::istringstream dbFile( DatabaseContents );
  FASTA::Reader< DNA> dbReader( dbFile );
  SequenceList< DNA > sequences;
  while( !dbReader.EndOfFile() ) {
    Sequence< DNA > seq;
    dbReader >> seq;
    sequences.push_back( std::move( seq ) );
  }

  db.Initialize( sequences );

  Sequence< DNA > query(
    "RF00807;mir-314;AAPT01020574.1/773332-773257   7222:Drosophila grimshawi",
    "UCGUAACUUGUGUGGCUUCGAAUGUACCUGGCUAAGGAAAGUUGGAUUUCCUAGGUAUUCGAGCCAAUAAGUUC"
    "GG" );

  SECTION( "Default" ) {
    GlobalSearch< DNA > gs( db, sp );
    auto hits = gs.Query( query );

    REQUIRE( hits.size() == 1 );
    REQUIRE( hits[ 0 ].target.identifier == "RF00807;mir-314;AFFE01007792.1/82767-82854   42026:Drosophila bipectinata" );
  }

  SECTION( "Min Identity" ) {
    sp.minIdentity = 0.9f;

    GlobalSearch< DNA > gs( db, sp );
    auto hits = gs.Query( query );

    REQUIRE( hits.size() == 0 );
  }

  SECTION( "Max Accepts" ) {
    sp.minIdentity = 0.6f;
    sp.maxAccepts = 2;

    GlobalSearch< DNA > gs( db, sp );
    auto hits = gs.Query( query );

    REQUIRE( hits.size() == 2 );
    std::vector< std::string > ids;
    for( auto &hit : hits ) {
      ids.push_back( hit.target.identifier );
    }
    REQUIRE( std::find( ids.begin(), ids.end(), "RF00807;mir-314;AAPU01011627.1/156896-156990   7230:Drosophila mojavensis" ) != ids.end() );
    REQUIRE( std::find( ids.begin(), ids.end(), "RF00807;mir-314;AFFE01007792.1/82767-82854   42026:Drosophila bipectinata" ) != ids.end() );
  }

  SECTION( "Strand support" ) {
    // our read goes in the "other" direction
    query = query.Reverse().Complement();

    SECTION( "Looking on plus strand (default)" ) {
      // no match, default is plus
      GlobalSearch< DNA > gs( db, sp );

      auto hits = gs.Query( query );
      REQUIRE( hits.size() == 0 );
    }

    SECTION( "Looking on minus strand" ) {
      // match, we are looking on the minus strand now
      sp.strand = DNA::Strand::Minus;
      GlobalSearch< DNA > gs( db, sp );

      auto hits = gs.Query( query );
      REQUIRE( hits.size() == 1 );
      REQUIRE( hits[ 0 ].strand == DNA::Strand::Minus );
    }

    SECTION( "Looking on both strands" ) {
      sp.strand = DNA::Strand::Both;
      GlobalSearch< DNA > gs( db, sp );

      auto hits = gs.Query( query );
      REQUIRE( hits.size() == 1 );
      REQUIRE( hits[ 0 ].strand == DNA::Strand::Minus );
    }
  }
}
