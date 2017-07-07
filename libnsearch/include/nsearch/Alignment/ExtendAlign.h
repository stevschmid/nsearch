#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#include "../Utils.h"

#define MAXINT INT_MAX/2 //prevent overflow
#define MININT -INT_MIN/2 //prevent underflow

typedef struct {
  int matchScore = 2;
  int mismatchScore = -4;

  int gapOpenScore = -20;
  int gapExtendScore = -2;

  int terminalGapOpenScore = -2;
  int terminalGapExtendScore = -1;

  inline int MatchScore( bool match ) const {
    return match ? matchScore : mismatchScore;
  };

  inline int GapOpenScore( bool terminal ) const {
    return terminal ? terminalGapOpenScore : gapOpenScore;
  }

  inline int GapExtendScore( bool terminal ) const {
    return terminal ? terminalGapExtendScore : gapExtendScore;
  }

  inline int GapScore( bool terminal, size_t length ) const {
    return GapOpenScore( terminal ) + length * GapExtendScore( terminal );
  }
} AlignmentParams;

enum class AlignmentDirection {
  forwards, backwards
};

class BandedAlign {
private:
  struct Cell {
    int score = MININT;

    struct {
      int score = MININT;
      int scoreToExtend = 0;
    } vGap;
  };
  using Cells = std::vector< Cell >;

  void PrintRow( size_t width ) {
    for( int i = 0; i < width; i++ ) {
      auto &c = mRow[ i ];
      if( c.score <= MININT ) {
        printf( "%5c", 'X' );
      }  else {
        printf( "%5d", c.score );
      }
    }
    printf("\n");
  }

  AlignmentParams mAP;
  Cells mRow;

public:
  int Align( const Sequence &A, const Sequence &B,
      size_t bandWidth,
      const AlignmentParams& ap = AlignmentParams(),
      AlignmentDirection dir = AlignmentDirection::forwards,
      size_t startA = 0, size_t startB = 0 )
  {
    int score;
    size_t x, y;
    size_t aIdx, bIdx;
    size_t width, height;

    // Calculate matrix width, depending on alignment
    // direction and length of sequences
    if( dir == AlignmentDirection::forwards ) {
      width = A.Length() - startA + 1;
      height = B.Length() - startB + 1;
    } else {
      width = startA + 1;
      height = startB + 1;
    }

    // Make sure we have enough cells
    if( mRow.capacity() < width ) {
      mRow = Cells( width * 1.5 );
    }

    // Initialize first row
    bool interiorX = startA > 0 && startA < A.Length();
    bool interiorY = startB > 0 && startB < B.Length();

    int bestScore = 0;
    mRow[ 0 ].score = 0;
    mRow[ 0 ].vGap.score = mAP.GapScore( !interiorY, 1 );
    mRow[ 0 ].vGap.scoreToExtend = mAP.GapExtendScore( !interiorY );

    for( x = 1; x < width && x <= bandWidth; x++ ) {
      mRow[ x ].score = mAP.GapScore( !interiorX, x );
      mRow[ x ].vGap.score = MININT;
    }
    PrintRow( width );

    size_t prevCenterX = 0;

    for( y = 1; y < height; y++ ) {

      int rowGap = MININT;
      int rowGapExtend = MININT;
      int score = MININT;

      size_t centerX = prevCenterX + 1;

      // Make sure we don't jump too far, so we can stil "connect the rows"
      size_t firstX = centerX > bandWidth ? ( centerX - bandWidth ) : 0;
      size_t lastX = std::min( centerX + bandWidth, width - 1 );

      if( y == height - 1 ) {
        lastX = width - 1;
      }

      int diagScore = MININT;
      if( firstX > 0 ) {
        diagScore = mRow[ firstX - 1 ].score;
        mRow[ firstX - 1 ].score = MININT;
      }

      for( x = firstX; x <= lastX; x++ ) {
        int colGap = mRow[ x ].vGap.score;
        int colGapExtend = mRow[ x ].vGap.scoreToExtend;

        aIdx = 0;
        bIdx = 0;
        if( x > 0 ) {
          aIdx = ( dir == AlignmentDirection::forwards ) ? startA + x - 1 : startA - x;
          bIdx = ( dir == AlignmentDirection::forwards ) ? startB + y - 1 : startB - y;
          // diagScore: score at col-1, row-1
          score = diagScore + ( DoNucleotidesMatch( A[ aIdx ], B[ bIdx ] ) ? mAP.matchScore : mAP.mismatchScore );
        }

        // select highest score
        //  - coming from diag (current),
        //  - coming from left (row)
        //  - coming from top (col)
        if( score < rowGap )
          score = rowGap;
        if( score < colGap )
          score = colGap;

        // Save the prev score at (x - 1) which
        // we will use to compute the diagonal score at (x)
        diagScore = mRow[ x ].score;

        // Record new score
        mRow[ x ].score = score;

        {
          // Horizontal gaps
          bool terminalGap = bIdx == 0 || bIdx + 1 == height - 1;
          int newGapOpen = terminalGap ? mAP.terminalGapOpenScore : mAP.gapOpenScore;
          int newGapExtend = terminalGap ? mAP.terminalGapExtendScore : mAP.gapExtendScore;

          int scoreNewGap = score + newGapOpen + newGapExtend;
          int scoreExtendingGap = rowGap + rowGapExtend;
          if( scoreNewGap > scoreExtendingGap ) {
            rowGap = scoreNewGap;
            rowGapExtend = newGapExtend;
          } else {
            rowGap = scoreExtendingGap;
          }
        }

        {
          // Vertical gaps
          bool terminalGap = aIdx == 0 || aIdx + 1 == width - 1;
          int newGapOpen = terminalGap ? mAP.terminalGapOpenScore : mAP.gapOpenScore;
          int newGapExtend = terminalGap ? mAP.terminalGapExtendScore : mAP.gapExtendScore;

          int scoreNewGap = score + newGapOpen + newGapExtend;
          int scoreExtendingGap = colGap + colGapExtend;
          if( scoreNewGap > scoreExtendingGap ) {
            mRow[ x ].vGap.score = scoreNewGap;
            mRow[ x ].vGap.scoreToExtend = newGapExtend;
          } else {
            mRow[ x ].vGap.score = scoreExtendingGap;
          }
        }

      }
      PrintRow( width );

      prevCenterX = centerX;
    }

    return bestScore;
  }
};

// Influenced by Blast's SemiGappedAlign function
class XDropExtendAlign {
private:
  struct Cell {
    int score = MININT;
    int scoreGap = MININT;
  };
  using Cells = std::vector< Cell >;

  void Print( const Cells& row ) {
    for( auto &c : row ) {
      if( c.score <= MININT ) {
        printf( "%5c", 'X' );
      }  else {
        printf( "%5d", c.score );
      }
    }
    printf("\n");
  }

  AlignmentParams mAP;
  Cells mRow;

public:
  XDropExtendAlign( const AlignmentParams& ap = AlignmentParams() )
    : mAP( ap )
  {
  }

  const AlignmentParams& AP() const {
    return mAP;
  }

  // Heavily influenced by Blast's SemiGappedAlign function
  int Extend( const Sequence &A, const Sequence &B,
      int xDrop,
      AlignmentDirection dir = AlignmentDirection::forwards,
      size_t startA = 0, size_t startB = 0,
      size_t *bestA = NULL, size_t *bestB = NULL )
  {
    int score;
    size_t x, y;
    size_t aIdx, bIdx;

    size_t width, height;

    if( dir == AlignmentDirection::forwards ) {
      width = A.Length() - startA + 1;
      height = B.Length() - startB + 1;
    } else {
      width = startA + 1;
      height = startB + 1;
    }

    if( mRow.capacity() < width ) {
      // Enlarge vector
      std::cout << "Enlarge row from " << mRow.capacity()
        << " to " << width << std::endl;
      mRow = Cells( width * 2 );
    }

    if( bestA ) {
      *bestA = startA;
    }
    if( bestB ) {
      *bestB = startB;
    }

    int bestScore = 0;
    mRow[ 0 ].score = 0; mRow[ 0 ].scoreGap = mAP.gapOpenScore + mAP.gapExtendScore;

    for( x = 1; x < width; x++ ) {
      score = mAP.gapOpenScore + x * mAP.gapExtendScore;

      if( score < -xDrop )
        break;

      mRow[ x ].score = score;
      mRow[ x ].scoreGap = MININT;
    }
    size_t rowSize = x;

    /* Print( mRow ); */

    size_t firstX = 0;

    for( y = 1; y < height; y++ ) {

      int rowGap = MININT;
      int score = MININT;
      int diagScore = MININT;

      size_t lastX = firstX;

      for( x = firstX; x < rowSize; x++ ) {
        int colGap = mRow[ x ].scoreGap;

        aIdx = 0;
        bIdx = 0;
        if( x > 0 ) {
          // diagScore: score at col-1, row-1

          if( dir == AlignmentDirection::forwards ) {
            aIdx = startA + x - 1;
            bIdx = startB + y - 1;
          } else {
            aIdx = startA - x;
            bIdx = startB - y;
          }

          /* printf( "x:%zu y:%zu %c == %c\n", x, y, A[ aIdx ], B[ bIdx ] ); */

          score = diagScore + ( DoNucleotidesMatch( A[ aIdx ], B[ bIdx ] ) ? mAP.matchScore : mAP.mismatchScore );
        }

        // select highest score
        //  - coming from diag (current),
        //  - coming from left (row)
        //  - coming from top (col)
        if( score < rowGap )
          score = rowGap;
        if( score < colGap )
          score = colGap;

        // mRow[ x ] right now points to the previous row, so use this
        // in the next iteration for the diagonal computation of (x, y )
        diagScore = mRow[ x ].score;

        if( bestScore - score > xDrop ) {
          // X-Drop test failed
          mRow[ x ].score = MININT;

          if( x == firstX ) {
            // Tighten left bound
            firstX++;
          }
        } else {
          lastX = x;

          // Check if we achieved new highscore
          if( score > bestScore ) {
            bestScore = score;
            if( bestA ) {
              *bestA = aIdx;
            }
            if( bestB ) {
              *bestB = bIdx;
            }
          }

          // Record new score
          mRow[ x ].score = score;
          mRow[ x ].scoreGap = std::max( score + mAP.gapOpenScore + mAP.gapExtendScore, colGap + mAP.gapExtendScore );
          rowGap = std::max( score + mAP.gapOpenScore + mAP.gapExtendScore, rowGap + mAP.gapExtendScore );
        }
      }

      if( firstX == rowSize ) {
        // All cells failed the X-Drop test
        // We are done 8)
        break;
      }

      if( lastX < rowSize - 1 ) {
        // Tighten right bound
        rowSize = lastX + 1;
      } else {
        // Extend row, since last checked column didn't fail X-Drop test
        while( rowGap >= ( bestScore - xDrop ) && rowSize < width ) {
          mRow[ rowSize ].score = rowGap;
          mRow[ rowSize ].scoreGap = rowGap + mAP.gapOpenScore + mAP.gapExtendScore;
          rowGap += mAP.gapExtendScore;
          rowSize++;
        }
      }

      /* Print( mRow ); */
    }

    return bestScore;
  };
};
