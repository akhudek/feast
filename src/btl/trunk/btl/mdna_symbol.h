/*
 *  mdna_symbol.h
 *  btl
 *
 *  Created by Alexander K. Hudek on 2008-08-06.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */
#ifndef BTL_MDNA_SYMBOL__
#define BTL_MDNA_SYMBOL__
#include <iostream>
#include <vector>
#include <stdexcept>
#include <inttypes.h>

namespace btl {
	class dna_symbol;
	class dna_gap_symbol;
	class dna_iupac_symbol;
	class mdna_iupac_symbol;
	class mdna_iupac_gap_symbol;

	class mdna_symbol {
		friend class btl::dna_symbol;
		friend class btl::dna_gap_symbol;
		friend class btl::dna_iupac_symbol;
		friend class btl::mdna_iupac_symbol;
		friend class btl::mdna_iupac_gap_symbol;

	protected:
		uint8_t value;

		inline uint8_t encode( char const c ) const {
			switch( c ) {
				case 'a':
					return 128;
				case 'A':
					return 0;
				case 't':
					return 129;
				case 'T':
					return 1;
				case 'c':
					return 130;
				case 'C':
					return 2;
				case 'g':
					return 131;
				case 'G':
					return 3;
				default:
					throw std::runtime_error("btl::mdna: unknown symbol");
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
				case 128:
					return 'a';
				case 129:
					return 't';
				case 130:
					return 'c';
				case 131:
					return 'g';
				default:
					throw std::runtime_error("btl::dna: unknown symbol");
			}
		}

	public:
		mdna_symbol();
		mdna_symbol(char c);
		mdna_symbol(mdna_symbol const &o);

		// conversion constructors
		mdna_symbol(dna_gap_symbol const &o);
		mdna_symbol(mdna_iupac_symbol const &o);

		// return an index to be used for an array
		inline  int index() const { return value&127; }

		// casting to an int is the same as asking for an index value
		inline operator int() const { return (int) value&127; }

		// dna sequences are never masked
		inline bool masked() const { return value&128; }

      // dna sequences are never masked
		inline void set_masked( bool m ) { value = (value&127) | (m << 7); }

		// expand ambiguity, for compatablitiy
		inline std::vector<mdna_symbol> expand_ambiguity() const { return std::vector<mdna_symbol>(1,*this); }

		bool operator==( mdna_symbol const &o ) const;

		bool operator!=( mdna_symbol const &o ) const;

		inline friend std::ostream &operator<<( std::ostream &out, mdna_symbol const &c );

		/** Return the compliment of the symbol. */
		inline mdna_symbol compliment() const {
			mdna_symbol r;
			r.value = ((value & 127) ^ 1) | (value & 128);
			return r;
		}
	};

	inline bool mdna_symbol::operator==( mdna_symbol const &o ) const { return (127&value) == (127&o.value); }

	inline bool mdna_symbol::operator!=( mdna_symbol const &o ) const { return (127&value) != (127&o.value); }

	inline  std::ostream &operator<<( std::ostream &out, btl::mdna_symbol const &c ) {
		out << c.decode(c.value);
		return out;
	}
}

#endif
