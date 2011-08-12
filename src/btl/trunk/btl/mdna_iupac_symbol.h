/*
 *  mdna_iupac_symbol.h
 *  btl
 *
 *  Created by Alexander K. Hudek on 2008-08-06.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */
#ifndef BTL_MDNA_IUPAC_SYMBOL__
#define BTL_MDNA_IUPAC_SYMBOL__
#include <iostream>
#include <vector>
#include <stdexcept>
#include <inttypes.h>

namespace btl {
	class mdna_symbol;
	class mdna_iupac_symbol;
	
	// representative sets
	extern std::vector< std::vector<mdna_iupac_symbol> > mdna_iupac_rsets;
	
	// here we construct the representative sets
	void mdna_iupac_construct_rsets();
	
	class mdna_iupac_symbol {
		friend class mdna_symbol;
		
	protected:
		uint8_t value;
		// high bit 8 is mask bit
		inline uint8_t encode( char const c ) const {
			switch( c ) {
				case 'g': return 128;
				case 'G': return 0;
				case 'c': return 129;
				case 'C': return 1;
				case 's': return 130;
				case 'S': return 2;
				case 't': return 131;
				case 'T': return 3;
				case 'k': return 132;
				case 'K': return 4;
				case 'y': return 133;
				case 'Y': return 5;
				case 'b': return 134;
				case 'B': return 6;
				case 'a': return 135;
				case 'A': return 7;
				case 'r': return 136;
				case 'R': return 8;
				case 'm': return 137;
				case 'M': return 9;
				case 'v': return 138;
				case 'V': return 10;
				case 'w': return 139;
				case 'W': return 11;
				case 'd': return 140;
				case 'D': return 12;
				case 'h': return 141;
				case 'H': return 13;
				case 'n': return 142;
				case 'N': return 14;
				default:
					throw std::runtime_error("btl::mdna_iupac_symbol: unknown symbol");
			}
		}
		
		inline char decode( uint8_t const b ) const {
			switch( b ) {
				case 0: return 'G';
				case 128: return 'g';
				case 1: return 'C';
				case 129: return 'c';
				case 2: return 'S';
				case 130: return 's';
				case 3: return 'T';
				case 131: return 't';
				case 4: return 'K';
				case 132: return 'k';
				case 5: return 'Y';
				case 133: return 'y';
				case 6: return 'B';
				case 134: return 'b';
				case 7: return 'A';
				case 135: return 'a';
				case 8: return 'R';
				case 136: return 'r';
				case 9: return 'M';
				case 137: return 'm';
				case 10: return 'V';
				case 138: return 'v';
				case 11: return 'W';
				case 139: return 'w';
				case 12: return 'D';
				case 140: return 'd';
				case 13: return 'H';
				case 141: return 'h';
				case 14: return 'N';
				case 142: return 'n';
				default:
					throw std::runtime_error("btl::mdna_iupac_symbol: unknown symbol");
			}
		}
		
	public:
		mdna_iupac_symbol();
		mdna_iupac_symbol(char c);
		mdna_iupac_symbol(mdna_iupac_symbol const &o);
		
		// conversion constructor
		mdna_iupac_symbol(mdna_symbol const &o);
		
		// return an index to be used for an array
		inline unsigned int index() const;
		
		// casting to an int is the same as asking for an index value
		inline operator int() const { return (int) value&127; }
		
		// is this symbol masked
		inline bool masked() const;
		
		// construct masked version of this symbol
		mdna_iupac_symbol mask() const;
		
		// expand ambiguity
		std::vector<mdna_iupac_symbol> const &represents() const;
		
		bool operator==( mdna_iupac_symbol const &o ) const;
		
		bool operator!=( mdna_iupac_symbol const &o ) const;
		
		friend mdna_iupac_symbol consensus( mdna_iupac_symbol const c1, mdna_iupac_symbol const c2 );
		
		friend std::ostream &operator<<( std::ostream &out, btl::mdna_iupac_symbol const &c );
		
	};
	
	inline std::vector<mdna_iupac_symbol> const &mdna_iupac_symbol::represents() const {
		if( mdna_iupac_rsets.empty() ) mdna_iupac_construct_rsets();
		return this->masked() ? mdna_iupac_rsets[this->index()+15] : mdna_iupac_rsets[this->index()];
	}
	
	// return an index to be used for an array
	inline unsigned int mdna_iupac_symbol::index() const { return value&127; }
	
	// dna sequences are never masked
	inline bool mdna_iupac_symbol::masked() const { return value&128; }

	// dna sequences are never masked
	inline mdna_iupac_symbol mdna_iupac_symbol::mask() const { 
		mdna_iupac_symbol newsymbol(*this);
		newsymbol.value |= 128;
		return newsymbol;
	}

	
	inline mdna_iupac_symbol consensus( mdna_iupac_symbol c1, mdna_iupac_symbol c2 ) {
		uint8_t mask = (c1.value&128) | (c2.value&128);
		c1.value &= 127;
		c2.value &= 127;
		c1.value++;
		c2.value++;
		mdna_iupac_symbol result;
		result.value = mask | ((c1.value & c2.value) == 0) ? (c1.value | c2.value) - 1 : (c1.value & c2.value) - 1;
		return result;
	}
	
	inline bool mdna_iupac_symbol::operator==( mdna_iupac_symbol const &o ) const { return (value&127) == (o.value&127); }
	
	inline bool mdna_iupac_symbol::operator!=( mdna_iupac_symbol const &o ) const { return (value&127) != (o.value&127); } 
	
	inline  std::ostream &operator<<( std::ostream &out, btl::mdna_iupac_symbol const &c ) {
		out << c.decode(c.value);
		return out;
	}
	
	
	
}
#endif

