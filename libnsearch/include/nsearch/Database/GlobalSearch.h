#pragma once

#include "BaseSearch.h"

#include "nsearch/Alignment/ExtendAlign.h"
#include "nsearch/Alignment/BandedAlign.h"

class GlobalSearch : public BaseSearch {
public:
  GlobalSearch( const Database &db, float minIdentity, int maxHits = 1, int maxRejects = 8 );
  QueryResult Query( const Sequence &query );

private:
  const Database &mDB;

  float mMinIdentity;
  int mMaxHits;
  int mMaxRejects;

  std::vector< size_t > mHits;
  ExtendAlign mExtendAlign;
  BandedAlign mBandedAlign;
};
