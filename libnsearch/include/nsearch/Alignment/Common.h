#pragma once

#include "../Sequence.h"
#include "../Utils.h"


#include <sstream>
#include <deque>

#define MAXINT INT_MAX/2 //prevent overflow
#define MININT -INT_MIN/2 //prevent underflow

enum class AlignmentDirection { fwd, rev };
