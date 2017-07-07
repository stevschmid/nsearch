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

  int terminalGapOpenScore = -2;
  int terminalGapExtendScore = -1;
} AlignmentParams;

enum class AlignExtendDirection { forwards, backwards };

class DP {
protected:
  struct Cell {
    int best = MININT;

    struct {
      int best = MININT;
      int extendingScore = 0;
    } vGap;
  };
  using Cells = std::vector< Cell >;

  AlignmentParams mAP;
  Cells mRow;
  AlignExtendDirection mDirection;
  size_t mWidth, mHeight, mStartA, mStartB;

  const Sequence *mA = NULL;
  const Sequence *mB = NULL;

  void PrintRow() {
    for( int i = 0; i < mWidth; i++ ) {
      auto &c = mRow[ i ];
      if( c.best <= MININT ) {
        printf( "%5c", 'X' );
      }  else {
        printf( "%5d", c.best );
      }
    }
    printf("\n");
  }

  inline int MatchScore( bool match ) const {
    return match ? mAP.matchScore : mAP.mismatchScore;
  };

  inline int GapOpenScore( bool terminal ) const {
    return terminal ? mAP.terminalGapOpenScore : mAP.gapOpenScore;
  }

  inline int GapExtendScore( bool terminal ) const {
    return terminal ? mAP.terminalGapExtendScore : mAP.gapExtendScore;
  }

  inline int GapScore( bool terminal, size_t length ) const {
    return GapOpenScore( terminal ) + length * GapExtendScore( terminal );
  }

  inline size_t MapXToSequenceA( size_t x ) const {
    return mDirection == AlignExtendDirection::forwards ? mStartA + x - 1 : mStartA - x;
  }

  inline size_t MapYToSequenceB( size_t y ) const {
    return mDirection == AlignExtendDirection::forwards ? mStartB + y - 1 : mStartB - y;
  }

  virtual void InitializeFirstRow() = 0;
  virtual bool ComputeRow( size_t y ) = 0;
  virtual int Traceback( size_t *outA = NULL, size_t *outB = NULL ) = 0;

public:
  DP( const AlignmentParams &ap = AlignmentParams() )
    : mAP( ap )
  {
  }

  int Align( const Sequence &A, const Sequence &B,
      AlignExtendDirection dir = AlignExtendDirection::forwards,
      size_t startA = 0, size_t startB = 0 )
  {
    mA = &A;
    mB = &B;

    mStartA = startA;
    mStartB = startB;

    mDirection = dir;

    if( dir == AlignExtendDirection::forwards ) {
      mWidth = A.Length() - startA + 1;
      mHeight = B.Length() - startB + 1;
    } else {
      mWidth = startA + 1;
      mHeight = startB + 1;
    }

    if( mRow.capacity() < mWidth ) {
      // Allocate new row if current one is insufficient in size
      mRow = Cells( mWidth * 1.5  );
    }

    InitializeFirstRow();
    for( size_t y = 1; y < mHeight; y++ ) {
      if( !ComputeRow( y ) )
        break;
    }
    return Traceback();
  }
};

class BandedDP : public DP {
protected:
  size_t mBandWidth;
  size_t mPreviousCenter;

  virtual void InitializeFirstRow() {
    size_t x;
    mRow[ 0 ].best = 0;
    mRow[ 0 ].vGap.best = GapScore( true, 1 );
    mRow[ 0 ].vGap.extendingScore = GapExtendScore( true );

    for( x = 1; x < mWidth && x <= mBandWidth; x++ ) {
      mRow[ x ].best = GapScore( true, x );
      mRow[ x ].vGap.best = MININT;
    }
    PrintRow();

    mPreviousCenter = 0;
  }

  virtual bool ComputeRow( size_t y ) {
    size_t aIdx, bIdx;

    int rowGap = MININT;
    int rowGapExtend = MININT;
    int score = MININT;

    size_t center = mPreviousCenter + 1;

    // Make sure we don't jump too far, so we can stil "connect the rows"
    size_t firstX = center > mBandWidth ? ( center - mBandWidth ) : 0;
    size_t lastX = std::min( center + mBandWidth, mWidth - 1 );

    if( y == mHeight - 1 ) {
      lastX = mWidth - 1;
    }

    int diagScore = MININT;
    if( firstX > 0 ) {
      diagScore = mRow[ firstX - 1 ].best;
      mRow[ firstX - 1 ].best = MININT;
    }

    for( size_t x = firstX; x <= lastX; x++ ) {
      int colGap = mRow[ x ].vGap.best;
      int colGapExtend = mRow[ x ].vGap.extendingScore;

      aIdx = 0;
      bIdx = 0;
      if( x > 0 ) {
        aIdx = MapXToSequenceA( x );
        bIdx = MapYToSequenceB( y );
        // diagScore: score at col-1, row-1
        bool match = DoNucleotidesMatch( (*mA)[ aIdx ], (*mB)[ bIdx ] );
        score = diagScore + MatchScore( match );
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
      diagScore = mRow[ x ].best;

      // Record new score
      mRow[ x ].best = score;

      {
        // Horizontal gaps
        bool terminal = bIdx == 0 || bIdx + 1 == mHeight - 1;

        int scoreNewGap = score + GapOpenScore( terminal ) + GapExtendScore( terminal );
        int scoreExtendingGap = rowGap + rowGapExtend;

        if( scoreNewGap > scoreExtendingGap ) {
          rowGap = scoreNewGap;
          rowGapExtend = GapExtendScore( terminal );
        } else {
          rowGap = scoreExtendingGap;
        }
      }

      {
        // Vertical gaps
        bool terminal = aIdx == 0 || aIdx + 1 == mWidth - 1;

        int scoreNewGap = score + GapOpenScore( terminal ) + GapExtendScore( terminal );
        int scoreExtendingGap = colGap + colGapExtend;

        if( scoreNewGap > scoreExtendingGap ) {
          mRow[ x ].vGap.best = scoreNewGap;
          mRow[ x ].vGap.extendingScore = GapExtendScore( terminal );
        } else {
          mRow[ x ].vGap.best = scoreExtendingGap;
        }
      }

    }
    PrintRow();
    mPreviousCenter = center;

    return true;
  }

  virtual int Traceback( size_t *outA = NULL, size_t *outB = NULL ) {
    return 0;
  }

public:
  BandedDP( size_t bandWidth,
      const AlignmentParams &ap = AlignmentParams() )
    : DP( ap ), mBandWidth( bandWidth )
  {
  }
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

    int bestScore = 0;
    mRow[ 0 ].score = 0;
    mRow[ 0 ].vGap.score = mAP.terminalGapOpenScore + mAP.terminalGapExtendScore;
    mRow[ 0 ].vGap.scoreToExtend = mAP.terminalGapExtendScore;

    for( x = 1; x < width && x <= bandWidth; x++ ) {
      score = mAP.terminalGapOpenScore + x * mAP.terminalGapExtendScore;
      mRow[ x ].score = score;
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
          aIdx = ( dir == AlignExtendDirection::forwards ) ? startA + x - 1 : startA - x;
          bIdx = ( dir == AlignExtendDirection::forwards ) ? startB + y - 1 : startB - y;
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
          bool terminalGap = bIdx == 0 || bIdx == height - 1;
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
          bool terminalGap = aIdx == 0 || aIdx == width - 1;
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
