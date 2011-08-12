#ifndef __BTL_ALGORITHM__
#define __BTL_ALGORITHM__

#include <algorithm>
#include <functional>
#include <vector>
#include <fstream>
#include <boost/lambda/lambda.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <btl/sequence.h>
#include <blitz/tinyvec.h>
#include <btl/gff_feature.h>

namespace btl {

   // Count letter frequencies in a sequence.
   template< typename ALPHA, typename CONTAINER >
   blitz::TinyVector<unsigned int, ALPHA::SIZE> symbol_frequencies( btl::sequence<ALPHA,CONTAINER> &seq ) {
      blitz::TinyVector<unsigned int, ALPHA::SIZE> count((unsigned int)0); 
      for( typename CONTAINER::iterator i = seq.data.begin(); i != seq.data.end(); i++ ) {
         count(i->index())++;
      }
      return count;
   }

   // Disambiguate an IUPAC character.
   template< typename URNG, typename IN_ALPHA, typename OUT_ALPHA = IN_ALPHA >
   struct disambiguate_iupac : public std::unary_function< typename IN_ALPHA::symbol, typename OUT_ALPHA::symbol > {
      private:
         // store a random number generator
         URNG &u01;

      public:
         // constructor
         disambiguate_iupac( URNG &rng ) : u01( rng ) {}
	
         // operator ()
         typename OUT_ALPHA::symbol operator()( typename IN_ALPHA::symbol const &b ) {
				boost::uniform_smallint<> pickone(0,b.represents().size()-1);
				return typename OUT_ALPHA::symbol( b.represents()[pickone(u01)] );
			}
   };

   // Test for gap.
   template< typename ALPHA >
   struct is_gap : public std::unary_function< typename ALPHA::symbol,  bool > {
      bool operator()( typename ALPHA::symbol const &b ) { return ALPHA::GAP == b; }
   };


	template< typename OUTPUT_ITERATOR >
	void read_gff_feature_set( std::istream &in, OUTPUT_ITERATOR out ) {
		while( 1 ) {
			std::string instr;
			while( instr.size() == 0 ) {
				if( in.eof() ) break;
				getline( in, instr );
				boost::algorithm::trim(instr);
				if( instr[0] == '#' ) instr.clear();
			}
			// nothing to load
			if( instr.size() == 0 ) break;
	
			gff_feature feature;
			feature.assign_from_string( instr );
			*out = feature;
			++out;
		}
	}

}

#endif

