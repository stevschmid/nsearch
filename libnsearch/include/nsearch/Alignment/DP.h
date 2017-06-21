#pragma once

#include <iostream>
#include <stdio.h>

#include <string>
#include <sstream>
#include <cassert>

#include "../Utils.h"

#define MININT -INT_MIN/2 //prevent underflow

class AlignmentParams {
public:
  int matchScore = 1;
  int mismatchScore = -2;

  int terminalGapOpenPenalty = 1;
  int terminalGapExtensionPenalty = 1;

  int interiorGapOpenPenalty = 10;
  int interiorGapExtensionPenalty = 1;

  int GapOpenScore( bool terminal ) const {
    return -( terminal ? terminalGapOpenPenalty : interiorGapOpenPenalty );
  }

  int GapExtensionScore( int length, bool terminal ) const {
    return length * -( terminal ? terminalGapExtensionPenalty : interiorGapExtensionPenalty );
  }

  int GapScore( int length, bool terminal ) const {
    return GapOpenScore( terminal ) + GapExtensionScore( length, terminal );
  }
};

class DPAlign {
public:
  DPAlign( const Sequence &A, const Sequence &B, AlignmentParams ap = AlignmentParams() )
    : mSequenceA( A ), mSequenceB( B ), mAP( ap )
  {
    // Todo: Intelligent cache
    mWidth = A.Length() + 1;
    mHeight = B.Length() + 1;
    mCells = new Cell[ mWidth * mHeight ];
  }

  virtual ~DPAlign() {
    delete[] mCells;
  }

  void DebugPrint( bool withArrows = false ) {
    for( size_t x = 0; x < mWidth; x++ ) {
      if( x == 0 ) {
        printf( "      " );
      } else {
        printf( "     %c", mSequenceA[ x - 1 ] );
      }
    }
    printf( "\n");

    for( size_t y = 0; y < mHeight; y++ ) {
      for( size_t x = 0; x < mWidth; x++ ) {
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

  virtual void ComputeMatrix() = 0;

  std::string Cigar() const {
    size_t x, y;
    if( !TracebackStartingPosition( x, y ) )
      return "";

    std::string aln;
    while( true ) {
      const Cell &cur = cell( x, y );
      if( cur.score <= MININT )
        break;

      if( cur.score == cur.hGap ) {
        // go left
        if( x == 0 )
          break;

        x--;
        aln += 'I';
      } else if( cur.score == cur.vGap ) {
        // go up
        if( y == 0 )
          break;

        aln += 'D';
        y--;
      } else {
        // go diagonal
        if( x == 0 || y == 0 )
          break;

        x--;
        y--;
        aln += 'M';
      }
    }

    std::stringstream cigar;
    char lastChar = 0;
    size_t lastCharCount = 0;
    for( auto it = aln.rbegin(); it != aln.rend(); it++ ) {
      if( *it != lastChar ) {
        if( lastChar ) {
          cigar << lastCharCount << lastChar;
        }
        lastChar = *it;
        lastCharCount = 0;
      }
      lastCharCount++;
    }
    if( lastChar ) {
      cigar << lastCharCount << lastChar;
    }

    return cigar.str();
  }

protected:
  virtual bool TracebackStartingPosition( size_t &x, size_t &y ) const = 0;

  void ComputeCell( size_t x, size_t y ) {
    Cell& cur = cell( x, y );

    if( x == 0 && y == 0 ) {
      cur.score = cur.vGap = cur.hGap = 0;
      return;
    }

    if( y == 0 ) {
      cur.score = cur.hGap = mAP.GapScore( x, true );
      return;
    }

    if( x == 0 ) {
      cur.score = cur.vGap = mAP.GapScore( y, true );
      return;
    }

    const Cell& leftCell = cell( x - 1, y );
    const Cell& upperCell = cell( x, y - 1 );
    const Cell& leftUpperCell = cell( x - 1, y - 1 );

    cur.hGap = std::max(
        leftCell.score + mAP.GapScore( 1, x == mWidth - 1 ),
        leftCell.hGap + mAP.GapExtensionScore( 1, x == mWidth - 1 ) );

    cur.vGap = std::max(
        upperCell.score + mAP.GapScore( 1, y == mHeight -1 ),
        upperCell.vGap + mAP.GapExtensionScore( 1, y == mHeight - 1 ) );

    int diagScore = leftUpperCell.score +
      ( DoNucleotidesMatch( mSequenceA[ x - 1 ], mSequenceB[ y - 1 ] ) ? mAP.matchScore : mAP.mismatchScore );
    cur.score = std::max( { diagScore, cur.hGap, cur.vGap } );
  }

  struct Cell {
    int score = MININT;
    int vGap = MININT;
    int hGap = MININT;
  };

  inline const Cell& cell( size_t x, size_t y ) const {
    assert( y >= 0 && y < mHeight );
    assert( x >= 0 && x < mWidth );
    return mCells[ y * mWidth + x ];
  }

  inline Cell& cell( size_t x, size_t y ) {
    assert( y >= 0 && y < mHeight );
    assert( x >= 0 && x < mWidth );
    return mCells[ y * mWidth + x ];
  }

  size_t mWidth, mHeight;
  Sequence mSequenceA, mSequenceB;
  Cell *mCells;
  AlignmentParams mAP;
};

class BandedGlobalAlign : public DPAlign
{
public:
  BandedGlobalAlign( const Sequence &A, const Sequence &B, AlignmentParams ap, size_t bandWidth )
    : DPAlign( A, B, ap ), mBandWidth( bandWidth )
  {
  }

  void ComputeMatrix() {
    for( size_t x = 0; x < mWidth; x++ ) {
      if( x > mBandWidth )
        break;

      ComputeCell( x, 0 );
    }

    size_t xFirst;
    size_t xLast;

    size_t lastCursor = 0;

    for( size_t y = 1; y < mHeight; y++ ) {

      size_t cursor = mWidth * float( y + 1 ) / float( mHeight );

      size_t xFirst = ( cursor > mBandWidth ) ? ( cursor - mBandWidth ) : 0;
      size_t xLast = ( cursor + mBandWidth < mWidth - 1 ) ? ( cursor + mBandWidth ) : ( mWidth - 1 );

      // Make sure we can connect to the previous row
      // (in case we jump far)
      if( xFirst > lastCursor ) {
        xFirst = lastCursor;
      }

      std::cout << "xfirst " << xFirst << " cursor " << cursor << " xlast " << xLast << std::endl;

      for( size_t x = xFirst; x <= xLast; x++ ) {
        ComputeCell( x, y );
      } // for x

      lastCursor = cursor;
    } // for y
  }

  bool TracebackStartingPosition( size_t &x, size_t &y ) const {
    x = mWidth - 1;
    y = mHeight - 1;
    return true;
  }

protected:
  size_t mBandWidth;
};
