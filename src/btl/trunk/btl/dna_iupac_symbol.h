/*
 *  dna_iupac_symbol.h
 *  btl
 *
 *  Created by Alexander K. Hudek on 2008-08-06.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */
#ifndef BTL_DNA_IUPAC_SYMBOL__
#define BTL_DNA_IUPAC_SYMBOL__
#include <iostream>
#include <vector>
#include <stdexcept>
#include <inttypes.h>

namespace btl {
	class mdna_symbol;
	
	class dna_iupac_symbol;
	
	// representative sets
	extern std::vector< std::vector<dna_iupac_symbol> > dna_iupac_rsets;

	// here we construct the representative sets
	void dna_iupac_construct_rsets();
	
	class dna_iupac_symbol {
	protected:
		uint8_t value;
		
		inline uint8_t encode( char const c ) const {
			switch( c ) {
				case 'g': case 'G': return 0;
				case 'c': case 'C': return 1;
				case 's': case 'S': return 2;
				case 't': case 'T': return 3;
				case 'k': case 'K': return 4;
				case 'y': case 'Y': return 5;
				case 'b': case 'B': return 6;
				case 'a': case 'A': return 7;
				case 'r': case 'R': return 8;
				case 'm': case 'M': return 9;
				case 'v': case 'V': return 10;
				case 'w': case 'W': return 11;
				case 'd': case 'D': return 12;
				case 'h': case 'H': return 13;
				case 'n': case 'N': return 14;
				default:
					throw std::runtime_error("btl::dna_iupac_symbol: unknown symbol");
			}
		}
		
		inline char decode( uint8_t const b ) const {
			switch( b ) {
				case 0: return 'G';
				case 1: return 'C';
				case 2: return 'S';
				case 3: return 'T';
				case 4: return 'K';
				case 5: return 'Y';
				case 6: return 'B';
				case 7: return 'A';
				case 8: return 'R';
				case 9: return 'M';
				case 10: return 'V';
				case 11: return 'W';
				case 12: return 'D';
				case 13: return 'H';
				case 14: return 'N';
					
				default:
					throw std::runtime_error("btl::dna_iupac_symbol: unknown symbol");
			}
		}
		
	public:
		dna_iupac_symbol();
		dna_iupac_symbol(char c);
		dna_iupac_symbol(dna_iupac_symbol const &o);
		
		// conversion constructor
		dna_iupac_symbol(mdna_symbol const &o);
		
		// return an index to be used for an array
		inline unsigned int index() const;
		
		// dna sequences are never masked
		inline bool masked() const;
		
		// expand ambiguity
		std::vector<dna_iupac_symbol> const &represents() const;
	
		bool operator==( dna_iupac_symbol const &o ) const;
		
		bool operator!=( dna_iupac_symbol const &o ) const;
		
		friend dna_iupac_symbol consensus( dna_iupac_symbol const c1, dna_iupac_symbol const c2 );
		
		friend std::ostream &operator<<( std::ostream &out, btl::dna_iupac_symbol const &c );
		
	};

	inline std::vector<dna_iupac_symbol> const &dna_iupac_symbol::represents() const {
		if( dna_iupac_rsets.empty() ) dna_iupac_construct_rsets();
		return dna_iupac_rsets[index()];
	}
			
	// return an index to be used for an array
	inline unsigned int dna_iupac_symbol::index() const { return value; }
		
	// dna sequences are never masked
	inline bool dna_iupac_symbol::masked() const { return false; }
		
	inline dna_iupac_symbol consensus( dna_iupac_symbol c1, dna_iupac_symbol c2 ) {
		c1.value++;
		c2.value++;
		dna_iupac_symbol result;
		result.value = ((c1.value & c2.value) == 0) ? (c1.value | c2.value) - 1 : (c1.value & c2.value) - 1;
		return result;
	}

	inline bool dna_iupac_symbol::operator==( dna_iupac_symbol const &o ) const { return value == o.value; }
	
	inline bool dna_iupac_symbol::operator!=( dna_iupac_symbol const &o ) const { return value != o.value; } 
	
	inline  std::ostream &operator<<( std::ostream &out, btl::dna_iupac_symbol const &c ) {
		out << c.decode(c.value);
		return out;
	}
	
	
	
}
#endif

