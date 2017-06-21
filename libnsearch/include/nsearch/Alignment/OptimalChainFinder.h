#pragma once

#include "Alignment.h"

#include <deque>
#include <map>

/*
 * Solving the two-dimensional chain problem
 * Adapted from 'Algorithms on Strings, Trees, and Sequences', Gusfield 1997
 * O(n log n)
 */
class OptimalChainFinder
{
  class Rect {
  public:
    size_t x1, x2;
    size_t y1, y2;
    int score;
    Rect *prev;

    Rect( size_t x1, size_t x2, size_t y1, size_t y2, size_t score )
      : x1( x1 ), x2( x2 ), y1( y1 ), y2( y2 ), score( score ), prev( NULL )
    {
    }
  };

  typedef std::shared_ptr< Rect > RectRef;

public:
  OptimalChainFinder( const SeedList &seeds ) {
    std::deque< RectRef > rects;
    std::multimap< size_t, RectRef > points;
    std::multimap< size_t, RectRef > solutions;

    auto greaterOrEqualThan = []( const std::multimap< size_t, RectRef >& set, size_t value ) {
      return set.lower_bound( value );
    };
    auto lessOrEqualThan = []( const std::multimap< size_t, RectRef >& set, size_t value ) {
      if( set.size() == 0 )
        return set.end();

      auto it = set.upper_bound( value );
      if( it == set.begin() )
        return set.end();

      return --it;
    };

    for( auto &seed : seeds ) {
      RectRef rect( new Rect( seed.s1, seed.s1 + seed.length, seed.s2, seed.s2 + seed.length, seed.length ) );
      rects.push_back( rect );
    }

    // Seed order should not matter:
    // For the same x-value, x2 (rectangle end) needs to come before x1 (rectangle start)
    // (since we need the previous rect in the x1 case)
    // Insertion order is guaranteed in C++11
    for( auto &rect : rects ) {
      points.insert( std::pair< size_t, RectRef >( rect->x2, rect ) );
    }
    for( auto &rect : rects ) {
      points.insert( std::pair< size_t, RectRef >( rect->x1, rect ) );
    }

    // Go through each point, from left to right
    for( auto &p : points ) {

      size_t px = p.first;
      RectRef rect = p.second;

      if( px == rect->x1 ) {
        // Left end of rectangle

        // Find the closest rectangle which can precede this one
        auto closest = lessOrEqualThan( solutions, rect->y1 );
        if( closest != solutions.end() ) {
          RectRef closestRect = (*closest).second;
          rect->prev = closestRect.get();
          rect->score += closestRect->score;
        }
      } else {
        // Right end of rectangle

        // Find competing, higher up rectangles
        auto competitor = greaterOrEqualThan( solutions, rect->y1 );
        if( competitor == solutions.end() || competitor->second->score < rect->score ) {
          solutions.insert( competitor, std::pair< size_t, RectRef >( rect->y2, rect ) );
        }

        // Delete competing solutions this one has beat (higher up but lower score)
        auto it = competitor;
        while( it != solutions.end() ) {
          if( (*it).second->score < rect->score ) {
            it = solutions.erase( it );
          } else {
            it++;
          }
        }
      }
    }

    auto bestSolution = std::max_element( solutions.begin(), solutions.end(), [](
          const std::pair< size_t, RectRef > &left,
          const std::pair< size_t, RectRef > &right )
    {
      return left.second->score < right.second->score;
    });

    if( bestSolution != solutions.end() ) {
      Rect *rect = bestSolution->second.get();
      mOptimalChain.clear();
      while( rect ) {
        mOptimalChain.push_back( Seed( rect->x1, rect->y1, rect->x2 - rect->x1 ) );
        rect = rect->prev;
      }
      std::reverse( mOptimalChain.begin(), mOptimalChain.end() );
    }
  };

  SeedList OptimalChain() const {
    return mOptimalChain;
  };

private:
  SeedList mOptimalChain;
};
