#pragma once

#include "BaseSearch.h"

#include "nsearch/Alignment/BandedAlign.h"
#include "nsearch/Alignment/ExtendAlign.h"

using Counter = unsigned short;

template < typename Alphabet >
class GlobalSearch : public BaseSearch< Alphabet > {
public:
  GlobalSearch( const Database< Alphabet >& db, const float minIdentity,
                const int maxHits = 1, const int maxRejects = 8 );
  HitList< Alphabet > Query( const Sequence< Alphabet >& query );

private:
  using BaseSearch< Alphabet >::mDB;

  float mMinIdentity;
  int   mMaxHits;
  int   mMaxRejects;

  std::vector< Counter >  mHits;
  ExtendAlign< Alphabet > mExtendAlign;
  BandedAlign< Alphabet > mBandedAlign;
};
