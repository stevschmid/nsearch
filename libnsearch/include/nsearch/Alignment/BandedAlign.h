#pragma once

#include "Cigar.h"
#include "Common.h"

#include <cassert>
#include <iostream>
#include <vector>

typedef struct BandedAlignParams {
  size_t bandwidth = 16;

  int matchScore    = 2;
  int mismatchScore = -4;

  int interiorGapOpenScore   = -20;
  int interiorGapExtendScore = -2;

  int terminalGapOpenScore   = -2;
  int terminalGapExtendScore = -1;
} BandedAlignParams;

class BandedAlign {
private:
  class Gap {
  private:
    int  mScore;
    bool mIsTerminal;

    const BandedAlignParams& mParams;

  public:
    Gap( const BandedAlignParams& params ) : mParams( params ) {
      Reset();
    }

    // Open new or extend existing
    void OpenOrExtend( int score, bool terminal, size_t length = 1 ) {
      int newGapScore = score;
      if( length > 0 ) {
        newGapScore += ( terminal ? mParams.terminalGapOpenScore
                                  : mParams.interiorGapOpenScore ) +
                       length * ( terminal ? mParams.terminalGapExtendScore
                                           : mParams.interiorGapExtendScore );
      }

      mScore += length * ( mIsTerminal ? mParams.terminalGapExtendScore
                                       : mParams.interiorGapExtendScore );

      if( newGapScore > mScore ) {
        mScore      = newGapScore;
        mIsTerminal = terminal;
      }
    }

    bool IsTerminal() const {
      return mIsTerminal;
    }

    int Score() const {
      return mScore;
    }

    void Reset() {
      mScore      = MININT;
      mIsTerminal = false;
    }
  };

  using Scores = std::vector< int >;
  using Gaps   = std::vector< Gap >;

  void PrintRow( size_t width ) {
    for( int i = 0; i < width; i++ ) {
      int score = mScores[ i ];
      if( score <= MININT ) {
        printf( "%5c", 'X' );
      } else {
        printf( "%5d", score );
      }
    }
    printf( "\n" );
  }

  Scores            mScores;
  Gaps              mVerticalGaps;
  CigarOps          mOperations;
  BandedAlignParams mParams;

public:
  BandedAlign( const BandedAlignParams& params = BandedAlignParams() )
      : mParams( params ) {}

  int Align( const Sequence& A, const Sequence& B, Cigar* cigar = NULL,
             AlignmentDirection dir = AlignmentDirection::fwd,
             size_t startA = 0, size_t startB = 0, size_t endA = -1,
             size_t endB = -1 ) {
    // Calculate matrix width, depending on alignment
    // direction and length of sequences
    // A will be on the X axis (width of matrix)
    // B will be on the Y axis (height of matrix)
    size_t width, height;

    size_t lenA = A.Length();
    size_t lenB = B.Length();

    if( endA == ( size_t ) -1 ) {
      endA = ( dir == AlignmentDirection::fwd ? lenA : 0 );
    }

    if( endB == ( size_t ) -1 ) {
      endB = ( dir == AlignmentDirection::fwd ? lenB : 0 );
    }

    if( startA > lenA )
      startA = lenA;
    if( startB > lenB )
      startB = lenB;
    if( endA > lenA )
      endA = lenA;
    if( endB > lenB )
      endB = lenB;

    width  = ( endA > startA ? endA - startA : startA - endA ) + 1;
    height = ( endB > startB ? endB - startB : startB - endB ) + 1;

    // Make sure we have enough cells
    if( mScores.capacity() < width ) {
      mScores = Scores( width * 1.5, MININT );
    }

    if( mVerticalGaps.capacity() < width ) {
      mVerticalGaps = Gaps( width * 1.5, mParams );
    }

    if( mOperations.capacity() < width * height ) {
      mOperations = CigarOps( width * height * 1.5 );
    }

    // Initialize first row
    size_t bw = mParams.bandwidth;

    bool fromBeginningA = ( startA == 0 || startA == lenA );
    bool fromBeginningB = ( startB == 0 || startB == lenB );

    bool fromEndA = ( endA == 0 || endA == lenA );
    bool fromEndB = ( endB == 0 || endB == lenB );

    mScores[ 0 ] = 0;

    mVerticalGaps[ 0 ].Reset();
    mVerticalGaps[ 0 ].OpenOrExtend( mScores[ 0 ], fromBeginningB );

    Gap horizontalGap( mParams );

    size_t x, y;
    for( x = 1; x < width; x++ ) {
      if( x > bw && height > 1 ) // only break on BW bound if B is not empty
        break;

      horizontalGap.OpenOrExtend( mScores[ x - 1 ], fromBeginningA );
      mScores[ x ]     = horizontalGap.Score();
      mOperations[ x ] = CigarOp::INSERTION;
      mVerticalGaps[ x ].Reset();
    }
    if( x < width ) {
      mScores[ x ] = MININT;
      mVerticalGaps[ x ].Reset();
    }
    /* PrintRow( width ); */

    // Row by row...
    size_t center = 1;
    bool   hitEnd = false;
    for( y = 1; y < height && !hitEnd; y++ ) {
      int score = MININT;

      // Calculate band bounds
      size_t leftBound =
        std::min( center > bw ? ( center - bw ) : 0, width - 1 );
      size_t rightBound = std::min( center + bw, width - 1 );

      /* // If we are in the last row, make sure we calculate up to the last
       * cell (for traceback) */
      /* if( y == height - 1 ) { */
      /*   rightBound = width - 1; */
      /* } */

      // Set diagonal score for first calculated cell in row
      int diagScore = MININT;
      if( leftBound > 0 ) {
        diagScore                = mScores[ leftBound - 1 ];
        mScores[ leftBound - 1 ] = MININT;
        mVerticalGaps[ leftBound - 1 ].Reset();
      }

      // Calculate row within the band bounds
      horizontalGap.Reset();
      for( x = leftBound; x <= rightBound; x++ ) {
        // Calculate diagonal score
        size_t aIdx = 0, bIdx = 0;
        bool   match;
        if( x > 0 ) {
          aIdx =
            ( dir == AlignmentDirection::fwd ) ? startA + x - 1 : startA - x;
          bIdx =
            ( dir == AlignmentDirection::fwd ) ? startB + y - 1 : startB - y;
          // diagScore: score at col-1, row-1
          match = DoNucleotidesMatch( A[ aIdx ], B[ bIdx ] );
          score =
            diagScore + ( match ? mParams.matchScore : mParams.mismatchScore );
        }

        // Select highest score
        //  - coming from diag (current),
        //  - coming from left (row)
        //  - coming from top (col)
        if( score < horizontalGap.Score() )
          score = horizontalGap.Score();

        Gap& verticalGap = mVerticalGaps[ x ];
        if( score < verticalGap.Score() )
          score = verticalGap.Score();

        // Save the prev score at (x - 1) which
        // we will use to compute the diagonal score at (x)
        diagScore = mScores[ x ];

        // Save new score
        mScores[ x ] = score;

        CigarOp op;
        if( score == horizontalGap.Score() ) {
          op = CigarOp::INSERTION;
        } else if( score == verticalGap.Score() ) {
          op = CigarOp::DELETION;
        } else {
          op = match ? CigarOp::MATCH : CigarOp::MISMATCH;
        }
        mOperations[ y * width + x ] = op;

        // Calculate potential gaps
        bool isTerminalA = ( x == 0 || x == width - 1 ) && fromEndA;
        bool isTerminalB = ( y == 0 || y == height - 1 ) && fromEndB;

        horizontalGap.OpenOrExtend( score, isTerminalB );
        verticalGap.OpenOrExtend( score, isTerminalA );
      }

      if( rightBound + 1 < width ) {
        mScores[ rightBound + 1 ] = MININT;
        mVerticalGaps[ rightBound + 1 ].Reset();
      }

      hitEnd = ( rightBound == leftBound );

      /* PrintRow( width ); */

      // Move one cell over for the next row
      center++;
    }

    // Backtrack
    if( cigar ) {
      size_t bx = x - 1;
      size_t by = y - 1;

      CigarEntry ce;
      cigar->Clear();
      while( bx != 0 || by != 0 ) {
        CigarOp op = mOperations[ by * width + bx ];
        cigar->Add( op );

        switch( op ) {
        case CigarOp::INSERTION:
          bx--;
          break;
        case CigarOp::DELETION:
          by--;
          break;
        case CigarOp::MATCH:
          bx--;
          by--;
          break;
        case CigarOp::MISMATCH:
          bx--;
          by--;
          break;
        default:
          assert( true );
          break;
        }
      }

      cigar->Reverse();
    }

    // Calculate score & cut corners
    int score = mScores[ x - 1 ];
    if( x == width ) {
      // We reached the end of A, emulate going down on B (vertical gaps)
      size_t remainingB  = height - y;
      Gap&   verticalGap = mVerticalGaps[ x - 1 ];
      verticalGap.OpenOrExtend( score, verticalGap.IsTerminal(), remainingB );
      score = verticalGap.Score();

      // Add tails to backtrack info
      if( cigar && remainingB ) {
        cigar->Add( CigarEntry( remainingB, CigarOp::DELETION ) );
      }
    } else if( y == height ) {
      // We reached the end of B, emulate going down on A (horizontal gaps)
      size_t remainingA = width - x;
      horizontalGap.OpenOrExtend( score, horizontalGap.IsTerminal(),
                                  remainingA );
      score = horizontalGap.Score();

      // Add tails to backtrack info
      if( cigar && remainingA ) {
        cigar->Add( CigarEntry( remainingA, CigarOp::INSERTION ) );
      }
    }

    if( cigar && dir == AlignmentDirection::rev ) {
      cigar->Reverse();
    }

    return score;
  }
};
