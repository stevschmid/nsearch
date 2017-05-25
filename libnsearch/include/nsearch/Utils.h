#pragma once

#include <string>
#include <ctype.h>

static std::string trim( const std::string &tstr ) {
  std::string str = tstr;
  str.erase( str.begin(), std::find_if_not(str.begin(), str.end(), [](char c){ return std::isspace(c); }) );
  str.erase( std::find_if_not(str.rbegin(), str.rend(), [](char c){ return std::isspace(c); }).base(), str.end() );
  return str;
}

static std::string toupper( const std::string &ustr ) {
  std::string str = ustr;
  for( auto &ch : str )
    ch = toupper( ch );
  return str;
}

