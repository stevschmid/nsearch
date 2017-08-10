#pragma once

#include <nsearch/Database/BaseSearch.h>

#include <string>

template < typename Alphabet >
extern bool Search( const std::string& queryPath,
                    const std::string& databasePath,
                    const std::string& outputPath,
                    const SearchParams< Alphabet >& searchParams );
