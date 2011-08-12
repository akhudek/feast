#ifndef BTL_DNA__
#define BTL_DNA__

#include <iostream>
#include <vector>
#include <stdexcept>
#include <btl/dna_symbol.h>
#include <boost/array.hpp>

namespace btl {
	class dna {
	public:
		// useful for template constructions
		static int const SIZE = 4;

		// reference our symbol type
		typedef dna_symbol symbol;

		// type for the alphabet storage
		typedef boost::array<dna_symbol,4> alphabet_type;

		// list of alphabet symbols
		static boost::array<dna_symbol,4> const alphabet;
		
		// here we encode constant versions of the symbols
		static dna_symbol const A, T, C, G;
   };
	
}

#endif
