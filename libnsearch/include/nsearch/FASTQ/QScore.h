#pragma once

#include <math.h>

namespace FASTQ
{
  static const int Q_MAX_SCORE = 41;
  static const int Q_MIN_ASCII_BASE = 33; // '!'
  static const int Q_MAX_ASCII_BASE = Q_MIN_ASCII_BASE + Q_MAX_SCORE; // 'J'

  // Calculate of posterior Q Scores as outlined by Edgar & Flyvbjerg (2015)
  class QScore
  {
  public:
    double ScoreToProbability( int q ) const {
      return pow( 10.0, -double( q ) / 10.0 );
    }

    int ProbabilityToScore( double p ) const {
      int q =  round( -10.0 * log10( p ) );
      if( q > Q_MAX_SCORE )
        q = Q_MAX_SCORE;
      return q;
    }

    int CalculatePosteriorScoreForMatch( int q1, int q2 ) const {
      return mPosteriorScoresForMatch[ q1 ][ q2 ];
    }

    int CalculatePosteriorScoreForMismatch( int q1, int q2 ) const {
      return mPosteriorScoresForMismatch[ q1 ][ q2 ];
    }

    // Singleton implementation
    static QScore &Instance() {
      // Safe in C++11
      static QScore instance;
      return instance;
    }
    QScore( QScore const& ) = delete;
    QScore( QScore&& ) = delete;
    QScore& operator=( QScore const& ) = delete;
    QScore& operator=( QScore && ) = delete;

  private:
    double mPosteriorScoresForMatch[ Q_MAX_SCORE + 1 ][ Q_MAX_SCORE + 1 ];
    double mPosteriorScoresForMismatch[ Q_MAX_SCORE + 1 ][ Q_MAX_SCORE + 1 ];

    void PrecomputeScores() {
      for( int qx = 0; qx <= Q_MAX_SCORE; qx++ ) {
        double px = ScoreToProbability( qx );

        for( int qy = 0; qy <= Q_MAX_SCORE; qy++ ) {
          double py = ScoreToProbability( qy );

          double pm = ( px * py / 3.0 )
            / ( 1.0 - px - py + 4.0 * px * py / 3.0 ); // (8) in paper

          double pmm;
          if( px < py ) {
            // px < py: merged base is X
            pmm = px * ( 1.0 - py / 3.0 )
              / ( px + py - 4.0 * px * py / 3.0 ); // (9) in paper
          } else {
            // py < px: merged base is Y
            pmm = py * ( 1.0 - px / 3.0 )
              / ( px + py - 4.0 * px * py / 3.0 ); // (9) in paper
          }

          mPosteriorScoresForMatch[ qx ][ qy ] = ProbabilityToScore( pm );
          mPosteriorScoresForMismatch[ qx ][ qy ] = ProbabilityToScore( pmm );
        }
      }
    }

    QScore() {
      PrecomputeScores();
    }
  };
}