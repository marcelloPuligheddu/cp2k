#include <vector>
#include <unordered_map>
#include <iostream>
#include "UniqueArray.h"
#include <functional>

using std::cout;
using std::endl;

bool linear_search( const std::vector<double>& data, double* values, int size, int* sp ){
   // I know it can be faster. If data is (indirectly) sorted it will be even faster to rewrite it
//   cout << " Looking for " << values[0] << " of size " << size << " among " << data.size() << endl;
   for( unsigned int i=0; i < data.size(); i++ ){
      bool found = false;
      bool is_possible = true;
      for( int s=0; s<size and is_possible ; s++ ){
         if ( values[s] != data[i+s] ){
            is_possible = false;
         }
      }
      if ( is_possible ){
//         cout << "Found " << i << endl;
         found = true;
         (*sp) = i;
         return found;
      }
   }
//   cout << " not found " << endl;
   return false;
}


__inline__ size_t hash_combine( size_t lhs, size_t rhs ) {
  return lhs ^ ( rhs + 0x517cc1b727220a95 + (lhs << 6) + (lhs >> 2) );
}

size_t hash_3doubles( const double* const x ){
   size_t hash = std::hash<double>{}(x[0]);
   hash = hash_combine( hash, std::hash<double>{}(x[1]) );
   hash = hash_combine( hash, std::hash<double>{}(x[2]) );
   return hash;
}

size_t hash_doubles( const double* const x, size_t l ){
   size_t hash_x = 0;
   for( size_t i=0; i < l; i++ ){
      size_t hash_i = std::hash<double>{}(x[i]);
      hash_x = hash_combine( hash_x, hash_i );
   }
   return hash_x;
}

UniqueArray::UniqueArray(){
   internal_buffer = std::vector<double>() ;
   hm = std::unordered_map<size_t,size_t>();
   count = 0;
}

// First horrible inefficient implementation
unsigned int UniqueArray::add( const double* const __restrict__ values, const int & size ){
/*
      if ( count % 1 == 0 ){
         cout << " Add iter # " << count ;
         for( int i=0; i < size; i++ ){ cout << " " << values[i] ; }
         cout << endl;
      }
      count++;
*/
      size_t h;
      if ( size == 3 ){
         h = hash_3doubles(values);
      } else {
         h = hash_doubles(values, size);
      }
      auto search = hm.find(h);
      if ( search != hm.end() ) { return search->second; }

      unsigned int prev_size = (unsigned int) internal_buffer.size();
      hm.insert({ h, prev_size });
      internal_buffer.insert(internal_buffer.end(), values, values+size );

      return prev_size;
/*
      unsigned int prev_size = (unsigned int) internal_buffer.size();
      int sp=0;
      bool found = linear_search( internal_buffer, values, size, &sp );
      if ( found ){ return sp ; }
      
      internal_buffer.insert(internal_buffer.end(), values, values+size );
      return prev_size;
*/
}

void UniqueArray::show_ia(){
   for ( unsigned int i=0; i < internal_buffer.size(); i++ ){
      std::cout << i << " " << internal_buffer[i] << std::endl ;
   }
}



