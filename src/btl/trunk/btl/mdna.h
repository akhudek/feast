#ifndef BTL_MDNA__
#define BTL_MDNA__
/*
 * This class represents the DNA alphabet but also carries additional masking information.
 * The lower and upper case symbols still represent the same base.
 *
 */
#include <iostream>
#include <vector>
#include <stdexcept>
#include <inttypes.h>
#include <boost/array.hpp>
#include <btl/mdna_symbol.h>

namespace btl {
   class mdna {
	public:
		// useful for template constructions
		static int const SIZE = 4;

		// reference our symbol type
		typedef mdna_symbol symbol;
		
		// type for the alphabet storage
		typedef boost::array<mdna_symbol,4> alphabet_type;

		// list of alphabet symbols
		static alphabet_type const alphabet;
		
		// here we encode constant versions of the symbols
		static mdna_symbol const A, T, C, G;
		
	};
		
}

#endif
