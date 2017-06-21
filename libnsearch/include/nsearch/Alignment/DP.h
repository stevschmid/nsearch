#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#define MAXINT INT_MAX/2 //prevent overflow
#define MININT -INT_MIN/2 //prevent underflow

typedef struct {
  int matchScore = 1;
  int mismatchScore = -2;

  int terminalGapOpenPenalty = 1;
  int terminalGapExtensionPenalty = 1;

  int interiorGapOpenPenalty = 10;
  int interiorGapExtensionPenalty = 1;

  int xDrop = 16;
} AlignmentParams;

class SemiGappedAlign
{
public:
  SemiGappedAlign( const std::string &A, const std::string &B, AlignmentParams ap = AlignmentParams() )
    : mWidth( A.length() + 1 ), mHeight( B.length() + 1 ),
      mSequenceA( A ), mSequenceB( B ),
      mAP( ap )
  {
    // Todo: Intelligent cache
    mCells = new Cell[ mWidth * mHeight ];
  }

  ~SemiGappedAlign() {
    delete[] mCells;
  }

  void Do() {
    size_t rowWidth;
    size_t x, y;

    cell( 0, 0 ).score = 0;

    auto gapOpenScore = [&]( bool terminal ) { return -( terminal ? mAP.terminalGapOpenPenalty : mAP.interiorGapOpenPenalty ); };
    auto gapExtensionScore = [&]( int length, bool terminal ) { return length * -( terminal ? mAP.terminalGapExtensionPenalty : mAP.interiorGapExtensionPenalty ); };
    auto gapScore = [&]( int length, bool terminal ) { return gapOpenScore( terminal ) + gapExtensionScore( length, terminal ); };

    for( x = 1; x < mWidth; x++ ) {
      Cell& cur = cell( x, 0 );
      int score = gapScore( x, true );
      if( score < -mAP.xDrop )
        break;
      cur.score = cur.hGap = score;
    }
    rowWidth = x;

    for( y = 1; y < mHeight; y++ ) {
      Cell& cur = cell( 0, y );
      int score = gapScore( y, true );
      if( score < -mAP.xDrop )
        break;
      cur.score = cur.vGap = score;
    }

    int bestScore = 0;

    size_t xStart = 1;
    size_t xEnd = mWidth;

    for( size_t y = 1; y < mHeight; y++ ) {

      for( size_t x = xStart; x < rowWidth; x++ ) {
        Cell& cur = cell( x, y );
        const Cell& leftCell = cell( x - 1, y );
        const Cell& upperCell = cell( x, y - 1 );
        const Cell& leftUpperCell = cell( x - 1, y - 1 );

        cur.hGap = std::max(
            leftCell.score + gapScore( 1, x == mWidth - 1 ),
            leftCell.hGap + gapExtensionScore( 1, x == mWidth - 1 ) );

        cur.vGap = std::max(
            upperCell.score + gapScore( 1, y == mHeight -1 ),
            upperCell.vGap + gapExtensionScore( 1, y == mHeight - 1 ) );

        int diagScore = leftUpperCell.score + ( mSequenceA[ x - 1 ] == mSequenceB[ y - 1 ] ? mAP.matchScore : mAP.mismatchScore );

        int score = std::max( { diagScore, cur.hGap, cur.vGap } );

        if( bestScore - score > mAP.xDrop ) {
          if( x == xStart ) {
            xStart++;
          }
        } else {
          xEnd = x;
          cur.score = score;
          bestScore = std::max( cur.score, bestScore );
        }
      }

      if( xEnd < rowWidth - 1 ) {
        // Shorter bounds next row
        rowWidth = xEnd + 1;
      } else {
        // This row is not finished yet, so extend it until x-drop test fails
        while( rowWidth < mWidth ) {
          Cell& cur = cell( rowWidth, y );
          const Cell& prev = cell( rowWidth - 1, y );
          int score = prev.hGap + gapExtensionScore( 1, x == mWidth - 1 );
          if( bestScore - score > mAP.xDrop ) {
            break;
          }

          cur.score = cur.hGap = score;
          rowWidth++;
        }
      }

    }
  }

  void DebugPrint( bool withArrows = false ) {
    for( int x = 0; x < mWidth; x++ ) {
      if( x == 0 ) {
        printf( "      " );
      } else {
        printf( "     %c", mSequenceA[ x - 1 ] );
      }
    }
    printf( "\n");

    for( int y = 0; y < mHeight; y++ ) {
      for( int x = 0; x < mWidth; x++ ) {
        if( x == 0 ) {
          char ch = ' ';
          if( y > 0 ) {
            ch = mSequenceB[ y - 1 ];
          }
          printf( "%c", ch );
        }

        if( cell( x, y ).score > MININT ) {
          printf( " %4d", cell( x, y ).score );

          if( withArrows ) {
            if( cell( x, y ).score == cell( x, y ).hGap ) {
              printf("←");
            } else if( cell( x, y ).score == cell( x, y ).vGap ) {
              printf("↑");
            } else {
              printf("↖");
            }
          } else  {
            printf(" ");
          }
        } else {
          printf( "    X " );
        }
      }

      printf( "\n" );
    }
  }

private:
  struct Cell {
    int score = MININT;
    int vGap = MININT;
    int hGap = MININT;
  };

  Cell& cell( size_t x, size_t y ) {
    assert( y >= 0 && y < mHeight );
    assert( x >= 0 && x < mWidth );
    return mCells[ y * mWidth + x ];
  }

  size_t mWidth, mHeight;
  std::string mSequenceA, mSequenceB;
  Cell *mCells;
  AlignmentParams mAP;
};
