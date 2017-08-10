#include <catch.hpp>

#include <nsearch/Alnout/Writer.h>
#include <nsearch/Alphabet/Protein.h>
#include <nsearch/Alphabet/DNA.h>

#include <sstream>
#include <cstring>

const char AlnoutOutputForProtein[] = R"(Query >query1
 %Id   TLen  Target
 75%      8  target50
 62%      8  target114

 Query 8aa >query1
Target 8aa >target50

Qry 1 LAFQGVRN 8
      :||||||.
Tgt 1 MAFQGVRS 8

8 cols, 6 ids (75.0%), 0 gaps (0.0%)

 Query 8aa >query1
Target 8aa >target114

Qry 1 LAFQGVRN 8
      || ||  |
Tgt 1 LAGQGSAN 8

8 cols, 5 ids (62.5%), 0 gaps (0.0%)

Query >query2
 %Id   TLen  Target
 91%     14  target1337

 Query 16aa >query2
Target 14aa >target1337

Qry  6 YFDEATGVCPF 16
       |||||||:|||
Tgt  1 YFDEATGICPF 11

11 cols, 10 ids (90.9%), 0 gaps (0.0%)

)";

const char AlnoutOutputForDNA[] = R"(Query >RF00966;mir-676;ABRQ01840532.1/340-428   9813:Procavia capensis (cape rock hyrax)
 %Id   TLen  Target
 85%     89  RF00966;mir-676;AAGV020395671.1/1356-1444   9361:Dasypus novemcinctus (nine-banded armadillo)

 Query 89nt >RF00966;mir-676;ABRQ01840532.1/340-428   9813:Procavia capensis (cape rock hyrax)
Target 89nt >RF00966;mir-676;AAGV020395671.1/1356-1444   9361:Dasypus novemcinctus (nine-banded armadillo)

Qry  1 CUUUGCCUGAACGCAAGACUCUUCAACCUCAGGACUUGCAGAAUUGGUAGAAUGCCGUCC 60
       | |  ||||||| || ||||||||||| |||||||||||||||||  | |||||||||||
Tgt  1 CGUCACCUGAACUCAUGACUCUUCAACUUCAGGACUUGCAGAAUUAAUGGAAUGCCGUCC 60

Qry 61 UAAGGUUGUUGAGUUCUGUGUUUGGAGGC 89
       |||||||||||||||||| ||||   |||
Tgt 61 UAAGGUUGUUGAGUUCUGCGUUUCUGGGC 89

89 cols, 76 ids (85.4%), 0 gaps (0.0%)

)";

TEST_CASE( "Alnout" ) {
  SECTION( "Protein" ) {
    auto entry1 = std::make_pair( Sequence< Protein >( "query1", "LAFQGVRN" ),
                                  HitList< Protein >( {
                                    { { "target50", "MAFQGVRS" }, "1X6M1X" },
                                    { { "target114", "LAGQGSAN" }, "4M3X1M" },
                                  } ) );
    auto entry2 =
      std::make_pair( Sequence< Protein >( "query2", "GGGGGYFDEATGVCPF" ),
                      HitList< Protein >( {
                        { { "target1337", "YFDEATGICPFQQQ" }, "5I7M1X3M3D" },
                      } ) );


    std::ostringstream oss;
    Alnout::Writer< Protein > writer( oss );

    writer << entry1;
    writer << entry2;

    REQUIRE( oss.str() == AlnoutOutputForProtein );
  }

  SECTION( "DNA" ) {
    auto entry = std::make_pair(
      Sequence< DNA >( "RF00966;mir-676;ABRQ01840532.1/340-428   9813:Procavia "
                       "capensis (cape rock hyrax)",
                       "CUUUGCCUGAACGCAAGACUCUUCAACCUCAGGACUUGCAGAAUUGGUAGA"
                       "AUGCCGUCCUAAGGUUGUUGAGUUCUGUGUUUGGAGGC" ),
      HitList< DNA >( {
        { { "RF00966;mir-676;AAGV020395671.1/1356-1444   9361:Dasypus "
            "novemcinctus (nine-banded armadillo)",
            "CGUCACCUGAACUCAUGACUCUUCAACUUCAGGACUUGCAGAAUUAAUGGAAUGCCGUCCUAAGGU"
            "UGUUGAGUUCUGCGUUUCUGGGC" },
          "1M1X1M2X7M1X2M1X11M1X17M2X1M1X29M1X4M3X3M" },
      } ) );

    std::ostringstream oss;
    Alnout::Writer< DNA > writer( oss );

    writer << entry;

    REQUIRE( oss.str() == AlnoutOutputForDNA );
  }
}
