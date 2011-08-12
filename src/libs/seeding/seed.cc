#include "seed.h"

//--------------------------------------------------------------------------
// Constructors
//--------------------------------------------------------------------------
seed::seed() {}

seed::seed( std::string const &s ) : p(s) {
	for( std::string::const_iterator i = p.begin(); i != p.end(); ++i )
		m.push_back( *i == '1' );
}

// functions
int seed::self_overlap() const {
   int max_overlap = 0;
   for( unsigned int i = 1; i < p.size(); i++ ) {
      // shift the seed
      int overlap = 0;
      for( unsigned int j = i, k = 0; j < p.size(); j++, k++ ) {
         if( (p[j] == '1') && (p[k] == '1') ) overlap++;
      } // for j,k
      if( overlap > max_overlap) max_overlap = overlap;

   } // for i
   return max_overlap;
}



