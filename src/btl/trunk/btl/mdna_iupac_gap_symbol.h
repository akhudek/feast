/*
 *  mdna_iupac_symbol.h
 *  btl
 *
 *  Created by Alexander K. Hudek on 2008-08-06.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */
#ifndef BTL_MDNA_IUPAC_GAP_SYMBOL__
#define BTL_MDNA_IUPAC_GAP_SYMBOL__
#include <iostream>
#include <vector>
#include <stdexcept>
#include <inttypes.h>

namespace btl {
	class mdna_symbol;
	class mdna_iupac_gap_symbol;
	
	// representative sets
	extern std::vector< std::vector<mdna_iupac_gap_symbol> > mdna_iupac_gap_rsets;
	
	// here we construct the representative sets
	void mdna_iupac_gap_construct_rsets();
	
	class mdna_iupac_gap_symbol {
		friend class mdna_symbol;
		
	protected:
		uint8_t value;
		// high bit 8 is mask bit
		inline uint8_t encode( char const c ) const {
			switch( c ) {
				case '-': return 0;
				case 'g': return 129;
				case 'G': return 1;
				case 'c': return 130;
				case 'C': return 2;
				case 's': return 131;
				case 'S': return 3;
				case 't': return 132;
				case 'T': return 4;
				case 'k': return 133;
				case 'K': return 5;
				case 'y': return 134;
				case 'Y': return 6;
				case 'b': return 135;
				case 'B': return 7;
				case 'a': return 136;
				case 'A': return 8;
				case 'r': return 137;
				case 'R': return 9;
				case 'm': return 138;
				case 'M': return 10;
				case 'v': return 139;
				case 'V': return 11;
				case 'w': return 140;
				case 'W': return 12;
				case 'd': return 141;
				case 'D': return 13;
				case 'h': return 142;
				case 'H': return 14;
				case 'n': return 143;
				case 'N': return 15;
				default:
					throw std::runtime_error("btl::mdna_iupac_symbol: unknown symbol");
			}
		}
		
		inline char decode( uint8_t const b ) const {
			switch( b ) {
				case 0: return '-';
				case 1: return 'G';
				case 129: return 'g';
				case 2: return 'C';
				case 130: return 'c';
				case 3: return 'S';
				case 131: return 's';
				case 4: return 'T';
				case 132: return 't';
				case 5: return 'K';
				case 133: return 'k';
				case 6: return 'Y';
				case 134: return 'y';
				case 7: return 'B';
				case 135: return 'b';
				case 8: return 'A';
				case 136: return 'a';
				case 9: return 'R';
				case 137: return 'r';
				case 10: return 'M';
				case 138: return 'm';
				case 11: return 'V';
				case 139: return 'v';
				case 12: return 'W';
				case 140: return 'w';
				case 13: return 'D';
				case 141: return 'd';
				case 14: return 'H';
				case 142: return 'h';
				case 15: return 'N';
				case 143: return 'n';
				default:
					throw std::runtime_error("btl::mdna_iupac_symbol: unknown symbol");
			}
		}
		
	public:
		mdna_iupac_gap_symbol();
		mdna_iupac_gap_symbol(char c);
		mdna_iupac_gap_symbol(mdna_iupac_gap_symbol const &o);
		
		// conversion constructor
		mdna_iupac_gap_symbol(mdna_symbol const &o);
		
		// return an index to be used for an array
		inline unsigned int index() const;
		
		// casting to an int is the same as asking for an index value
		inline operator int() const { return (int) value&127; }
		
		// is this symbol masked
		inline bool masked() const;
		
		// construct masked version of this symbol
		mdna_iupac_gap_symbol mask() const;
		
		// expand ambiguity
		std::vector<mdna_iupac_gap_symbol> const &represents() const;
		
		bool operator==( mdna_iupac_gap_symbol const &o ) const;
		
		bool operator!=( mdna_iupac_gap_symbol const &o ) const;
		
		friend mdna_iupac_gap_symbol consensus( mdna_iupac_gap_symbol const c1, mdna_iupac_gap_symbol const c2 );
		
		friend std::ostream &operator<<( std::ostream &out, btl::mdna_iupac_gap_symbol const &c );
		
	};
	
	inline std::vector<mdna_iupac_gap_symbol> const &mdna_iupac_gap_symbol::represents() const {
		if( mdna_iupac_gap_rsets.empty() ) mdna_iupac_gap_construct_rsets();
		return this->masked() ? mdna_iupac_gap_rsets[this->index()+16] : mdna_iupac_gap_rsets[this->index()];
	}
	
	// return an index to be used for an array
	inline unsigned int mdna_iupac_gap_symbol::index() const { return value&127; }
	
	// dna sequences are never masked
	inline bool mdna_iupac_gap_symbol::masked() const { return value&128; }
	
	// dna sequences are never masked
	inline mdna_iupac_gap_symbol mdna_iupac_gap_symbol::mask() const { 
		mdna_iupac_gap_symbol newsymbol(*this);
		newsymbol.value |= 128;
		return newsymbol;
	}
	
	
	inline mdna_iupac_gap_symbol consensus( mdna_iupac_gap_symbol c1, mdna_iupac_gap_symbol c2 ) {
		uint8_t mask = (c1.value&128) | (c2.value&128);
		c1.value &= 127;
		c2.value &= 127;
		mdna_iupac_gap_symbol result;
		result.value = mask | ((c1.value & c2.value) == 0) ? c1.value | c2.value : c1.value & c2.value;
		return result;
	}
	
	inline bool mdna_iupac_gap_symbol::operator==( mdna_iupac_gap_symbol const &o ) const { return (value&127) == (o.value&127); }
	
	inline bool mdna_iupac_gap_symbol::operator!=( mdna_iupac_gap_symbol const &o ) const { return (value&127) != (o.value&127); } 
	
	inline  std::ostream &operator<<( std::ostream &out, btl::mdna_iupac_gap_symbol const &c ) {
		out << c.decode(c.value);
		return out;
	}
	
	
	
}
#endif

