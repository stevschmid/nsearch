#pragma once

#include "../Utils.h"

#define MAXINT INT_MAX/2 //prevent overflow
#define MININT -INT_MIN/2 //prevent underflow

enum class AlignmentDirection {
  forwards, backwards
};

typedef struct {
  int score;

  size_t a1, a2;
  size_t b1, b2;

  std::string ops;
} Alignment;

class BandedAlignParams
{
public:
  size_t bandwidth = 16;

  int matchScore = 2;
  int mismatchScore = -4;

  int interiorGapOpenScore = -20;
  int interiorGapExtendScore = -2;

  int terminalGapOpenScore = -2;
  int terminalGapExtendScore = -1;
};

class BandedAlign {
private:
  class Gap {
  private:
    int mScore;
    bool mIsTerminal;

    int mTerminalGapScore, mInteriorGapScore;
    int mTerminalGapExtendScore, mInteriorGapExtendScore;

  public:
    Gap( const BandedAlignParams &params )
      :
        mTerminalGapScore( params.terminalGapOpenScore + params.terminalGapExtendScore ),
        mTerminalGapExtendScore( params.terminalGapExtendScore ),

        mInteriorGapScore( params.interiorGapOpenScore + params.interiorGapExtendScore ),
        mInteriorGapExtendScore( params.interiorGapExtendScore )
    {
      Reset();
    }

    void OpenOrExtend( int score, bool terminal ) {
      int newGapScore = score + ( terminal ? mTerminalGapScore : mInteriorGapScore );
      mScore += ( mIsTerminal ? mTerminalGapExtendScore : mInteriorGapExtendScore );

      if( newGapScore > mScore ) {
        mScore = newGapScore;
        mIsTerminal = terminal;
      }
    }

    int Score() {
      return mScore;
    }

    void Reset() {
      mScore = MININT;
      mIsTerminal = false;
    }
  };

  using Scores = std::vector< int >;
  using Gaps = std::vector< Gap >;
  using Operations = std::vector< char >;

  void PrintRow( size_t width ) {
    for( int i = 0; i < width; i++ ) {
      int score = mScores[ i ];
      if( score <= MININT ) {
        printf( "%5c", 'X' );
      }  else {
        printf( "%5d", score );
      }
    }
    printf("\n");
  }

  Scores mScores;
  Gaps mVerticalGaps;
  Operations mOperations;
  BandedAlignParams mParams;

public:
  BandedAlign( const BandedAlignParams& params = BandedAlignParams() )
    : mParams( params )
  {

  }

  Alignment Align( const Sequence &A, const Sequence &B,
      bool backtrack = false,
      size_t startA = 0, size_t startB = 0,
      AlignmentDirection dir = AlignmentDirection::forwards )
  {
    Alignment aln;

    // Calculate matrix width, depending on alignment
    // direction and length of sequences
    // A will be on the X axis (width of matrix)
    // B will be on the Y axis (height of matrix)
    size_t width, height;

    size_t lenA = A.Length();
    size_t lenB = B.Length();

    if( dir == AlignmentDirection::forwards ) {
      width = lenA - startA + 1;
      height = lenB - startB + 1;
    } else {
      width = startA + 1;
      height = startB + 1;
    }

    // Make sure we have enough cells
    if( mScores.capacity() < width ) {
      mScores = Scores( width * 1.5, MININT );
    }

    if( mVerticalGaps.capacity() < width ) {
      mVerticalGaps = Gaps( width * 1.5, mParams );
    }

    if( mOperations.capacity() < width * height ) {
      mOperations = Operations( width * height * 1.5 );
    }

    // Initialize first row
    size_t bw = mParams.bandwidth;

    bool fromBeginningA = ( startA == 0 || startA == lenA );
    bool fromBeginningB = ( startB == 0 || startB == lenB );

    mScores[ 0 ] = 0;
    mVerticalGaps[ 0 ].OpenOrExtend( mScores[ 0 ], fromBeginningB );

    Gap horizontalGap( mParams );
    for( size_t x = 1; x < width && x <= bw; x++ ) {
      horizontalGap.OpenOrExtend( mScores[ x - 1 ], fromBeginningA );
      mScores[ x ] = horizontalGap.Score();
      mVerticalGaps[ x ].Reset();
      mOperations[ x ] = 'D';
    }
    PrintRow( width );

    // Row by row...
    size_t center = 0;
    for( size_t y = 1; y < height; y++ ) {
      int score = MININT;

      // Calculate band bounds
      size_t leftBound = std::min( center > bw ? ( center - bw ) : 0, width - 1 );
      size_t rightBound = std::min( center + bw, width - 1 );

      // If we are in the last row, make sure we calculate up to the last cell (for traceback)
      if( y == height - 1 ) {
        rightBound = width - 1;
      }

      // Set diagonal score for first calculated cell in row
      int diagScore = MININT;
      if( leftBound > 0 ) {
        diagScore = mScores[ leftBound - 1 ];
        mScores[ leftBound - 1 ] = MININT;
      }

      // Calculate row within the band bounds
      horizontalGap.Reset();
      for( size_t x = leftBound; x <= rightBound; x++ ) {
        // Calculate diagonal score
        size_t aIdx = 0, bIdx = 0;
        if( x > 0 ) {
          aIdx = ( dir == AlignmentDirection::forwards ) ? startA + x - 1 : startA - x;
          bIdx = ( dir == AlignmentDirection::forwards ) ? startB + y - 1 : startB - y;
          // diagScore: score at col-1, row-1
          bool match = DoNucleotidesMatch( A[ aIdx ], B[ bIdx ] );
          score = diagScore + ( match ? mParams.matchScore : mParams.mismatchScore );
        }

        // Select highest score
        //  - coming from diag (current),
        //  - coming from left (row)
        //  - coming from top (col)
        if( score < horizontalGap.Score() )
          score = horizontalGap.Score();

        Gap &verticalGap = mVerticalGaps[ x ];
        if( score < verticalGap.Score() )
          score = verticalGap.Score();

        // Save the prev score at (x - 1) which
        // we will use to compute the diagonal score at (x)
        diagScore = mScores[ x ];

        // Save new score
        mScores[ x ] = score;

        char op;
        if( score == horizontalGap.Score() ) {
          op = 'D';
        } else if( score == verticalGap.Score() ) {
          op = 'I';
        } else {
          op = 'M';
        }
        mOperations[ y * width + x ] = op;

        // Calculate potential gaps
        bool isTerminalA = ( x == 0 || x == width - 1 );
        bool isTerminalB = ( y == 0 || y == height - 1 );

        horizontalGap.OpenOrExtend( score, isTerminalB );
        /* if( y == 1 ){ */
        /*   std::cout << score << std::endl; */
        /*   std::cout << isTerminalB << std::endl; */
        /*   std::cout << x << " " << horizontalGap.Score() << std::endl; */
        /* } */
        verticalGap.OpenOrExtend( score, isTerminalA );
      }
      PrintRow( width );

      // Move one cell over for the next row
      center++;
    }

    // Backtrack
    if( backtrack ) {
      size_t x = width - 1;
      size_t y = height - 1;
      while( x != 0 || y != 0 ) {
        char op = mOperations[ y * width + x ];
        aln.ops.push_back( op );
        switch( op ) {
          case 'D': x--; break;
          case 'I': y--; break;
          case 'M': x--; y--; break;
        }
      }
      std::reverse( aln.ops.begin(), aln.ops.end() );
    }

    return aln;
  }
};
