#pragma once

enum class AlignmentDirection {
  forwards, backwards
};

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

  inline size_t Bandwidth() const {
    return bandwidth;
  }

  inline int MatchScore( bool match ) const {
    return match ? matchScore : mismatchScore;
  };

  inline int GapOpenScore( bool terminal ) const {
    return terminal ? terminalGapOpenScore : interiorGapOpenScore;
  }

  inline int GapExtendScore( bool terminal ) const {
    return terminal ? terminalGapExtendScore : interiorGapExtendScore;
  }

  inline int GapScore( bool terminal, size_t length ) const {
    return GapOpenScore( terminal ) + length * GapExtendScore( terminal );
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

  Cells mRow;
  BandedAlignParams mParams;

public:
  int Align( const Sequence &A, const Sequence &B,
      const BandedAlignParams& params = BandedAlignParams(),
      size_t startA = 0, size_t startB = 0,
      AlignmentDirection dir = AlignmentDirection::forwards
  {
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
    if( mRow.capacity() < width ) {
      mRow = Cells( width * 1.5 );
    }

    // Initialize first row
    size_t bw = params.BandWidth();

    bool terminalX = startA == 0 || startA == lenA;
    bool terminalY = startB == 0 || startB == lenB;

    int bestScore = 0;
    mRow[ 0 ].score = 0;
    mRow[ 0 ].vGap.score = params.GapScore( terminalY, 1 );
    mRow[ 0 ].vGap.scoreToExtend = params.GapExtendScore( terminalY );

    for( size_t x = 1; x < width && x <= bw; x++ ) {
      mRow[ x ].score = params.GapScore( terminalX, x );
      mRow[ x ].vGap.score = MININT;
    }

    // Row by row...
    size_t center = 0;
    for( y = 1; y < height; y++ ) {
      int score = MININT;

      int hGapScore = MININT;
      int hGapExtend = MININT;

      // Calculate band bounds
      size_t leftBound = center > bw ? ( center - bw ) : 0;
      size_t rightBound = std::min( center + bw, width - 1 );

      // If we are in the last row, make sure we calculate up to the last cell (for traceback)
      if( y == height - 1 ) {
        leftBound = width - 1;
      }

      // Set diagonal score for first calculated cell in row
      int diagScore = MININT;
      if( leftBound ) {
        diagScore = mRow[ leftBound - 1 ].score;
        mRow[ leftBound - 1 ].score = MININT;
      }

      // Calculate row within the band bounds
      aIdx = 0;
      bIdx = 0;

      for( x = leftBound; x <= rightBound; x++ ) {
        int vGapScore = mRow[ x ].vGap.score;
        int vGapExtend = mRow[ x ].vGap.scoreToExtend;

        // Calculate diagonal score
        if( x > 0 ) {
          aIdx = ( dir == AlignmentDirection::forwards ) ? startA + x - 1 : startA - x;
          bIdx = ( dir == AlignmentDirection::forwards ) ? startB + y - 1 : startB - y;
          // diagScore: score at col-1, row-1
          bool match = DoNucleotidesMatch( A[ aIdx ], B[ bIdx ] );
          score = diagScore + mParams.MatchScore( match );
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

        // Save new score
        mRow[ x ].score = score;

        // Calculate potential gaps
        {
          // Horizontal gaps

          // Score when extending the current (hypothetical) gap
          int scoreExtendingGap = rowGap + rowGapExtend;

          // Horizontal gaps
          bool terminal = bIdx == 0 || bIdx == lenB - 1;
          int scoreNewGap = score + mParams.GapCost( terminal, 1 );

          if( scoreNewGap > scoreExtendingGap ) {
            rowGap = scoreNewGap;
            rowGapExtend = newGapExtend;
          } else {
            rowGap = scoreExtendingGap;
          }
        }

        // TODO: Continue...

        {
          // Vertical gaps
          bool terminal = aIdx == 0 || aIdx + 1 == width - 1;
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

      // Move one cell over for the next row
      center++;
    }

    return bestScore;
  }
};
