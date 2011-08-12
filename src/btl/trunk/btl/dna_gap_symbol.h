/*
 *  dna_gap_symbol.h
 *  btl
 *
 *  Created by Alexander K. Hudek on 2008-08-06.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */
#ifndef BTL_DNA_GAP_SYMBOL__
#define BTL_DNA_GAP_SYMBOL__
#include <iostream>
#include <vector>
#include <stdexcept>
#include <inttypes.h>

namespace btl {
	class mdna_symbol;
	class dna_symbol;
	
	class dna_gap_symbol {
		friend class mdna_symbol;
		
	protected:
		uint8_t value;
		
		inline uint8_t encode( char const c ) const {
			switch( c ) {
				case 'a':
				case 'A':
					return 0;
				case 't':
				case 'T':
					return 1;
				case 'c':
				case 'C':
					return 2;
				case 'g':
				case 'G':
					return 3;
				case '-':
					return 4;
				default:
					throw std::runtime_error("btl::dna_gap_symbol: unknown symbol");
			}
		}
		
		inline char decode( uint8_t const b ) const {
			switch( b ) {
				case 0:
					return 'A';
				case 1:
					return 'T';
				case 2:
					return 'C';
				case 3:
					return 'G';
				case 4:
					return '-';
				default:
					throw std::runtime_error("btl::dna_gap_symbol: unknown symbol");
			}
		}
		
	public:
		dna_gap_symbol();
		dna_gap_symbol(char c);
		dna_gap_symbol(dna_gap_symbol const &o);
		dna_gap_symbol(btl::dna_symbol const &o);
		
		// conversion constructor
		dna_gap_symbol(mdna_symbol const &o);
		
		// return an index to be used for an array
		inline int index() const { return value; }
		
		// casting to an int is the same as asking for an index value
		inline operator int() const { return (int) value; }
		
		// dna sequences are never masked
		inline bool masked() const { return false; }
		
		// expand ambiguity, for compatablitiy
		inline std::vector<dna_gap_symbol> expand_ambiguity() const { return std::vector<dna_gap_symbol>(1,*this); }
		
		bool operator==( dna_gap_symbol const &o ) const;
		
		bool operator!=( dna_gap_symbol const &o ) const;
		
		friend std::ostream &operator<<( std::ostream &out, btl::dna_gap_symbol const &c );
		
	};
	
	inline bool dna_gap_symbol::operator==( dna_gap_symbol const &o ) const { return value == o.value; }
	
	inline bool dna_gap_symbol::operator!=( dna_gap_symbol const &o ) const { return value != o.value; } 
	
	inline  std::ostream &operator<<( std::ostream &out, btl::dna_gap_symbol const &c ) {
		out << c.decode(c.value);
		return out;
	}
	
}
#endif

