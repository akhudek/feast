#ifndef BTL_DNA_IUPAC__
#define BTL_DNA_IUPAC__

#include <inttypes.h>
#include <stdexcept>
#include <boost/array.hpp>
#include <btl/dna_iupac_symbol.h>

namespace btl {
	
	class dna_iupac {
	private:
		static int const SIZE = 15;
		
	public:
		// reference our symbol type
		typedef dna_iupac_symbol symbol;
		
		// list of alphabet symbols
		static boost::array<dna_iupac_symbol,15> const alphabet;
		
		// here we encode constant versions of the symbols
		static dna_iupac_symbol const G, C, S, T, K, Y, B, A, R, M, V, W, D, H, N;
		
   };
	
}

#endif
