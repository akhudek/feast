#ifndef BTL_MDNA_IUPAC_GAP__
#define BTL_MDNA_IUPAC_GAP__

#include <inttypes.h>
#include <stdexcept>
#include <boost/array.hpp>
#include <btl/mdna_iupac_gap_symbol.h>

namespace btl {
	
	class mdna_iupac_gap {
	private:
		static int const SIZE = 16;
		
	public:
		// reference our symbol type
		typedef mdna_iupac_gap_symbol symbol;
		
		// alphabet type
		typedef boost::array<mdna_iupac_gap_symbol,16> const alphabet_type;
		
		// list of alphabet symbols
		static alphabet_type const alphabet;
		
		// here we encode constant versions of the symbols
		static mdna_iupac_gap_symbol const G, C, S, T, K, Y, B, A, R, M, V, W, D, H, N, GAP;
		
   };
	
}

#endif
