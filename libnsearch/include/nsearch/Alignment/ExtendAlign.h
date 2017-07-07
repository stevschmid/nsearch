#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#include "../Utils.h"

#define MAXINT INT_MAX/2 //prevent overflow
#define MININT -INT_MIN/2 //prevent underflow

typedef struct { int matchScore = 2;
  int mismatchScore = -4;

  int gapOpenScore = -20;
  int gapExtendScore = -2;
} AlignmentParams;

enum class AlignExtendDirection { forwards, backwards };

class BandedAlign {
private:
  struct Cell {
    int score = MININT;
    int scoreGap = MININT;

    struct {
      int score = MININT;
      int scoreToExtend = 0;
    } vGap;
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
  BandedAlign( const AlignmentParams& ap = AlignmentParams() )
    : mAP( ap )
  {
  }

  int Align( const Sequence &A, const Sequence &B,
      size_t bandWidth,
      AlignExtendDirection dir = AlignExtendDirection::forwards,
      size_t startA = 0, size_t startB = 0 )
  {
    int score;
    size_t x, y;
    size_t aIdx, bIdx;

    size_t width, height;

    if( dir == AlignExtendDirection::forwards ) {
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
      mRow = Cells( width  );
    }

    auto SeqPosA = [&]( size_t x ) {
      return dir == AlignExtendDirection::forwards ? startA + x - 1 : startA - x;
    };
    auto SeqPosB = [&]( size_t y ) {
      return dir == AlignExtendDirection::forwards ? startB + y - 1 : startB - y;
    };

    int bestScore = 0;
    mRow[ 0 ].score = 0;
    /* mRow[ 0 ].scoreGap = mAP.gapOpenScore + mAP.gapExtendScore; */
    mRow[ 0 ].vGap.score = mAP.gapOpenScore + mAP.gapExtendScore;
    mRow[ 0 ].vGap.scoreToExtend = mAP.gapExtendScore;

    for( x = 1; x < width && x <= bandWidth; x++ ) {
      score = mAP.gapOpenScore + x * mAP.gapExtendScore;
      mRow[ x ].score = score;
      /* mRow[ x ].scoreGap = MININT; */
      mRow[ x ].vGap.score = MININT;
    }
    Print( mRow );
    size_t rowSize = x;

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
          aIdx = SeqPosA( x );
          bIdx = SeqPosB( y );
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

        int scoreNewVGap = score + mAP.gapOpenScore + mAP.gapExtendScore;
        int scoreExtendVGap = colGap + colGapExtend;

        if( scoreNewVGap > scoreExtendVGap ) {
          mRow[ x ].vGap.score = scoreNewVGap;
          mRow[ x ].vGap.scoreToExtend = mAP.gapExtendScore;
        } else {
          mRow[ x ].vGap.score = scoreExtendVGap;
        }

        int scoreNewHGap = score + mAP.gapOpenScore + mAP.gapExtendScore;
        int scoreExtendHGap = rowGap + rowGapExtend;
        if( scoreNewHGap > scoreExtendHGap ) {
          rowGap = scoreNewHGap;
          rowGapExtend = mAP.gapExtendScore;
        } else {
          rowGap = scoreExtendHGap;
        }
      }
      Print( mRow );

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
      AlignExtendDirection dir = AlignExtendDirection::forwards,
      size_t startA = 0, size_t startB = 0,
      size_t *bestA = NULL, size_t *bestB = NULL )
  {
    int score;
    size_t x, y;
    size_t aIdx, bIdx;

    size_t width, height;

    if( dir == AlignExtendDirection::forwards ) {
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

          if( dir == AlignExtendDirection::forwards ) {
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
