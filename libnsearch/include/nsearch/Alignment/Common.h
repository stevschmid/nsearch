#pragma once

#include "../Utils.h"

#define MAXINT INT_MAX/2 //prevent overflow
#define MININT -INT_MIN/2 //prevent underflow

enum class AlignmentDirection {
  forwards, backwards
};

using Cigar = std::string;
