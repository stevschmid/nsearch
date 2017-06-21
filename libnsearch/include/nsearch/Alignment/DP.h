#pragma once

#include <iostream>
#include <stdio.h>

#include <set>
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

  int Score() const {
    size_t x, y;
    if( !TracebackStartingPosition( x, y ) )
      return MININT;

    return cell( x, y ).score;
  }

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

class GuidedBandedGlobalAlign : public DPAlign
{
  class GuidePoint {
  public:
    size_t x, y;

    bool operator<( const GuidePoint &other ) const {
      if( y < other.y )
        return true;
      if( y > other.y )
        return false;
      return x < other.x;
    }

    GuidePoint( size_t x = 0, size_t y = 0 )
      : x( x ), y ( y )
    {
    }
  };

public:
  GuidedBandedGlobalAlign( const Sequence &A, const Sequence &B, AlignmentParams ap, size_t bandWidth, const SeedList &chain = SeedList() )
    : DPAlign( A, B, ap ), mBandWidth( bandWidth )
  {
    for( auto &pos : chain ) {
      mGuidePoints.insert( GuidePoint( pos.s1, pos.s2 ) );
      mGuidePoints.insert( GuidePoint( pos.s1 + pos.length, pos.s2 + pos.length ) );
    }

    if( !mGuidePoints.empty() ) {
      // Start position must intersect with matrix
      auto &fgp = *mGuidePoints.begin();
      if( fgp.x != 0 && fgp.y != 0 ){
        size_t offset = std::min( fgp.x, fgp.y );
        mGuidePoints.insert( GuidePoint( fgp.x - offset, fgp.y - offset ) );
      }

      // End position must intersect with matrix
      auto &lgp = *mGuidePoints.rbegin();
      if( lgp.x != mWidth - 1 && lgp.y != mHeight - 1 ) {
        size_t offset = std::min( mWidth - 1 - lgp.x, mHeight - 1 - lgp.y );
        mGuidePoints.insert( GuidePoint( lgp.x + offset, lgp.y + offset ) );
      }
    } else {
      // Straight diagonal
      mGuidePoints.insert( GuidePoint( 0, 0 ) );
      mGuidePoints.insert( GuidePoint( mWidth - 1, mHeight - 1 ) );
    }
  }

  void ComputeMatrix() {
    if( mGuidePoints.empty() )
      return;

    size_t yStart = (*mGuidePoints.begin()).y;
    size_t yEnd = (*mGuidePoints.rbegin()).y;

    size_t lastCursor = 0;

    for( size_t y = yStart; y <= yEnd; y++ ) {
      // Current guide point
      auto currentIt = mGuidePoints.upper_bound( GuidePoint( y, 0 ) );
      if( currentIt != mGuidePoints.begin() ) {
        currentIt--;
      }
      const GuidePoint& currentGP = *currentIt;

      // Next guide point
      auto nextIt = mGuidePoints.upper_bound( currentGP );
      assert( nextIt != mGuidePoints.end() );
      const GuidePoint& nextGP = *nextIt;

      float ratio = float( y - currentGP.y ) / float( nextGP.y - currentGP.y );
      size_t cursor = currentGP.x + ratio * ( nextGP.x - currentGP.x);

      size_t xFirst = ( cursor > mBandWidth ) ? ( cursor - mBandWidth ) : 0;
      size_t xLast = ( cursor + mBandWidth < mWidth - 1 ) ? ( cursor + mBandWidth ) : ( mWidth - 1 );

      // Make sure we can connect to the previous row
      // (in case we jump far)
      if( xFirst > lastCursor )
        xFirst = lastCursor;
      lastCursor = cursor;

      /* std::cout << "Y " << y */
      /*   << " CurrentGP " << "(" << currentGP.x << "," << currentGP.y << ")" */
      /*   << " NextGP " << "(" << nextGP.x << "," << nextGP.y << ")" */
      /*   << std::endl; */
      /* std::cout << "Ratio " << ratio << std::endl; */
      /* std::cout << "Cursor " << cursor << std::endl; */
      for( size_t x = xFirst; x <= xLast; x++ ) {
        ComputeCell( x, y );
      }
    } // for y
  }

  bool TracebackStartingPosition( size_t &x, size_t &y ) const {
    // Needleman Wunsch
    x = mWidth - 1;
    y = mHeight - 1;
    return true;
  }

private:
  size_t mBandWidth;
  std::multiset< GuidePoint > mGuidePoints;
};
