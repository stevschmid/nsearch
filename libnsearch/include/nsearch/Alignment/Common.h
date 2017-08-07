#pragma once

#include "../Sequence.h"
#include "../Utils.h"

#include <deque>
#include <sstream>

#define MAXINT INT_MAX / 2  // prevent overflow
#define MININT -INT_MIN / 2 // prevent underflow

enum class AlignmentDirection { Forward, Reverse };
