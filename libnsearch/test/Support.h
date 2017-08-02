#pragma once

#include <nsearch/Database/Kmers.h>
#include <string>

static Kmer Kmerify( const std::string& a ) {
  Kmer kmer = 0;
  size_t counter = 0;
  for( auto &ch : a ) {
    char val;
    switch( ch ) {
      case 'A': val = 0b00; break;
      case 'C': val = 0b01; break;
      case 'G': val = 0b11; break;
      default: val = 0b10; /// T, U
    }
    kmer |= val << counter;
    counter += 2;
  }
  return kmer;
}

