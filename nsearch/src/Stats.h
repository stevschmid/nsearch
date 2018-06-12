#pragma once

#include <atomic>
#include <chrono>

class Stats {
public:
  std::atomic< size_t > numProcessed;
  std::atomic< size_t > numMerged;
  std::atomic< size_t > mergedReadsTotalLength;

  Stats() : numProcessed( 0 ), numMerged( 0 ), mergedReadsTotalLength( 0 ) {}

  double MeanMergedLength() const {
    return float( mergedReadsTotalLength ) / numMerged;
  }

  void StartTimer() {
    mTimerStart = std::chrono::steady_clock::now();
  }

  void StopTimer() {
    mTimerStop = std::chrono::steady_clock::now();
  }

  double ElapsedMillis() const {
    auto diff = mTimerStop - mTimerStart;
    auto elapsedMs =
      std::chrono::duration_cast< std::chrono::milliseconds >( diff ).count();
    return elapsedMs;
  }

private:
  std::chrono::time_point< std::chrono::steady_clock > mTimerStart, mTimerStop;
};

extern Stats gStats;
