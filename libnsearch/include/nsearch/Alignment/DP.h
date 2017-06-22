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

  int terminalGapOpenPenalty = 0;
  int terminalGapExtensionPenalty = 0;

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

class Alignment {
public:
  Alignment()
    : mScore( 0 )
  {
  }

  Alignment( const Sequence &A, const Sequence &B, int score, const std::string &cigar )
    : mSequenceA( A ), mSequenceB( B ), mScore( score ), mCigar( cigar )
  {
  }

  int Score() const {
    return mScore;
  }

  std::string Cigar() const {
    return mCigar;
  }

private:
  Sequence mSequenceA, mSequenceB;
  int mScore;
  std::string mCigar;
};

class DPAlign {
public:
  DPAlign( AlignmentParams ap = AlignmentParams() )
    : mAP( ap ), mWidth( 0 ), mHeight( 0 ), mCells( NULL ), mNumCellsReserved( 0 )
  {
  }

  int Align( const Sequence &A, const Sequence &B, Alignment *aln = NULL ) {
    mSequenceA = A;
    mSequenceB = B;

    // Todo: Intelligent cache
    mWidth = A.Length() + 1;
    mHeight = B.Length() + 1;

    // Enlarge cells if we need more
    // (Reserving memory is a heavy task)
    if( mNumCellsReserved < mWidth * mHeight ) {
      std::cout << "Resize cells " << mNumCellsReserved << " -> " << (mWidth * mHeight) << std::endl;

      if( mCells )
        delete[] mCells;

      mNumCellsReserved = mWidth * mHeight;
      mCells = new Cell[ mNumCellsReserved ];
    }

    // Reset matrix
    Cell *cell = mCells;
    while( cell < mCells + mNumCellsReserved ) {
      cell->score = cell->vGap = cell->hGap = MININT;
      cell++;
    }

    // Now Compute Matrix
    ComputeMatrix();

    // Populate aln
    int score = Score();

    if( aln ) {
      *aln = Alignment( A, B, score, Cigar() );
    }

    return score;
  }

  virtual ~DPAlign() {
    if( mCells ) {
      delete[] mCells;
    }
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

    bool terminalRow = ( y == 0 || y == mHeight - 1 );
    bool terminalCol = ( x == 0 || x == mWidth - 1 );

    cur.hGap = std::max(
        leftCell.score + mAP.GapScore( 1, terminalRow ),
        leftCell.hGap + mAP.GapExtensionScore( 1, terminalRow ) );

    cur.vGap = std::max(
        upperCell.score + mAP.GapScore( 1, terminalCol ),
        upperCell.vGap + mAP.GapExtensionScore( 1, terminalCol ) );

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
  size_t mNumCellsReserved;
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
  GuidedBandedGlobalAlign( size_t bandWidth, AlignmentParams ap = AlignmentParams() )
    : DPAlign( ap ), mBandWidth( bandWidth )
  {
  }

  int AlignAlongChain( const Sequence &A, const Sequence &B, const SeedList &chain = SeedList(), Alignment *aln = NULL ) {
    mGuidePoints.clear();

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
      if( lgp.x != A.Length() && lgp.y != B.Length() ) {
        size_t offset = std::min( A.Length() - lgp.x, B.Length() - lgp.y );
        mGuidePoints.insert( GuidePoint( lgp.x + offset, lgp.y + offset ) );
      }
    } else {
      // Straight diagonal
      mGuidePoints.insert( GuidePoint( 0, 0 ) );
      mGuidePoints.insert( GuidePoint( A.Length(), B.Length() ) );
    }

    return Align( A, B, aln );
  }

  void ComputeMatrix() {
    if( mGuidePoints.empty() )
      return;

    const GuidePoint& firstGP = *mGuidePoints.begin();
    const GuidePoint& lastGP = *mGuidePoints.rbegin();

    size_t lastCursor = 0;

    // Calculate from 0, 0 to intersection point (first GP)
    if( firstGP.x == 0 ) {
      // Intersects at Y-axis
      for( size_t y = 0; y < firstGP.y; y++ ) {
        ComputeCell( 0, y );
      }
    } else {
      // Intersects at X-axis
      for( size_t x = 0; x < firstGP.x; x++ ) {
        ComputeCell( x, 0 );
      }
    }

    for( size_t y = firstGP.y; y <= lastGP.y; y++ ) {
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

    // Calculate from last GP to end of matrix
    if( lastGP.x == mWidth - 1 ) {
      // Intersects at Y-axis
      for( size_t y = lastGP.y; y < mHeight; y++ ) {
        ComputeCell( mWidth - 1, y );
      }
    } else {
      // Intersects at X-axis
      for( size_t x = lastGP.x; x < mWidth; x++ ) {
        ComputeCell( x, mHeight - 1 );
      }
    }
  }

  bool TracebackStartingPosition( size_t &tx, size_t &ty ) const {
    /* size_t bestX = 0, bestY = 0; */
    /* int bestScore = MININT; */

    /* for( size_t x = 0; x < mWidth; x++ ) { */
    /*   const Cell &cur = cell( x, mHeight - 1 ); */
    /*   if( cur.score > bestScore ) { */
    /*     bestScore = cur.score; */
    /*     bestX = x; */
    /*     bestY = mHeight - 1; */
    /*   } */
    /* } */

    /* for( size_t y = 0; y < mHeight; y++ ) { */
    /*   const Cell &cur = cell( mWidth - 1, y ); */
    /*   if( cur.score > bestScore ) { */
    /*     bestScore = cur.score; */
    /*     bestX = mWidth - 1; */
    /*     bestY = y; */
    /*   } */
    /* } */

    /* tx = bestX; */
    /* ty = bestY; */

    /* return true; */

    // Needleman-Wunsch
    tx = mWidth - 1;
    ty = mHeight - 1;

    return true;
  }

private:
  size_t mBandWidth;
  std::multiset< GuidePoint > mGuidePoints;
};
