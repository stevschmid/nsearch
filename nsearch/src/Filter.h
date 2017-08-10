#pragma once

#include <string>

extern bool DoFilter( const std::string& inputPath,
                      const std::string& outputPath,
                      const float        maxExpectedErrors );
