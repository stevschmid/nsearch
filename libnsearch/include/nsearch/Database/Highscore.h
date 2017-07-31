#pragma once

class Highscore {
  class Entry {
  public:
    size_t id    = 0;
    size_t score = 0;

    bool operator<( const Entry& other ) const {
      return score < other.score;
    }
  };

public:
  Highscore( const size_t numHighestEntriesToKeep ) : mLowestScore( 0 ) {
    mEntries.resize( numHighestEntriesToKeep );
  }

  // score is assumed to increase for every id
  void Set( const size_t id, const size_t score ) {
    if( score < mLowestScore )
      return;

    auto it = std::find_if(
      mEntries.begin(), mEntries.end(),
      [id]( const Entry& candidate ) { return id == candidate.id; } );

    if( it == mEntries.end() ) {
      it = std::find_if(
        mEntries.begin(), mEntries.end(),
        [score]( const Entry& candidate ) { return score > candidate.score; } );
    }

    if( it != mEntries.end() ) {
      it->id    = id;
      it->score = score;

      mLowestScore =
        std::min_element( mEntries.begin(), mEntries.end() )->score;
    }
  }

  std::vector< Entry > EntriesFromTopToBottom() const {
    std::vector< Entry > topToBottom = mEntries;
    std::sort( topToBottom.begin(), topToBottom.end(),
               []( const Entry& a, const Entry& b ) { return !( a < b ); } );
    return topToBottom;
  }

private:
  size_t               mLowestScore;
  std::vector< Entry > mEntries;
};
