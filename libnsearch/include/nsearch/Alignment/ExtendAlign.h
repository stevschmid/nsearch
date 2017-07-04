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

  int xDrop = MAXINT;
} AlignmentParams;

// Influenced by Blast's SemiGappedAlign function
class ExtendAlign {
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
  enum class ExtendDirection { forwards, backwards };

  ExtendAlign( const AlignmentParams& ap = AlignmentParams() )
    : mAP( ap )
  {
  }

  // Heavily influenced by Blast's SemiGappedAlign function
  int Extend( const Sequence &A, const Sequence &B,
      ExtendDirection dir = ExtendDirection::forwards,
      size_t startA = 0, size_t startB = 0,
      size_t *bestA = NULL, size_t *bestB = NULL )
  {
    int score;
    size_t x, y;
    size_t aIdx, bIdx;

    size_t width, height;

    if( dir == ExtendDirection::forwards ) {
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
    mRow[ 0 ].score = 0;
    mRow[ 0 ].scoreGap = mAP.gapOpenScore + mAP.gapExtendScore;

    for( x = 1; x < width; x++ ) {
      score = mAP.gapOpenScore + x * mAP.gapExtendScore;

      if( score < -mAP.xDrop )
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

          if( dir == ExtendDirection::forwards ) {
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

        // Save the prev score at (x - 1) which
        // we will use to compute the diagonal score at (x)
        diagScore = mRow[ x ].score;

        if( bestScore - score > mAP.xDrop ) {
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
        while( rowGap >= ( bestScore - mAP.xDrop ) && rowSize < width ) {
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
