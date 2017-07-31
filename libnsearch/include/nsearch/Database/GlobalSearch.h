#pragma once

#include "BaseSearch.h"

#include "nsearch/Alignment/BandedAlign.h"
#include "nsearch/Alignment/ExtendAlign.h"

using Counter = unsigned short;

class GlobalSearch : public BaseSearch {
public:
  GlobalSearch( const Database& db, const float minIdentity,
                const int maxHits = 1, const int maxRejects = 8 );
  HitList Query( const Sequence& query );

private:
  const Database& mDB;

  float mMinIdentity;
  int   mMaxHits;
  int   mMaxRejects;

  std::vector< Counter > mHits;
  ExtendAlign            mExtendAlign;
  BandedAlign            mBandedAlign;
};
