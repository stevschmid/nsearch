#pragma once

#include <string>

extern bool Search( const std::string& queryPath,
                    const std::string& databasePath,
                    const std::string& outputPath, float minIdentity,
                    int maxAccepts, int maxRejects );
