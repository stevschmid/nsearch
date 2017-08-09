#pragma once

#include <string>

template < typename Alphabet >
extern bool Search( const std::string& queryPath,
                    const std::string& databasePath,
                    const std::string& outputPath, const float minIdentity,
                    const int maxAccepts, const int maxRejects );
