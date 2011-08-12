#ifndef TRIPLET_TABLE_H__
#define TRIPLET_TABLE_H__
#include <blitz/tinyvec.h>
#include "dna_sequence.h"

template<typename DATA = unsigned int>
class triplet_table {
private:
blitz::TinyVector< DATA,64> data;

public:
triplet_table() { data = (DATA)0; }

inline DATA &operator()( dna_alpha::symbol x2, dna_alpha::symbol x1, dna_alpha::symbol x ) {
	return data(16*x2.index() + 4*x1.index() + x.index());
}

inline triplet_table<DATA> operator+( triplet_table<DATA> const &other ) {
	triplet_table<DATA> newtri = *this;
	newtri.data += other.data;
	return newtri;
}

inline triplet_table<DATA> &operator+=( triplet_table<DATA> const &other ) {
	data += other.data;
	return *this;
}

inline triplet_table<DATA> operator-( triplet_table<DATA> const &other ) {
	triplet_table<DATA> newtri = *this;
	newtri.data -= other.data;
	return newtri;
}

friend std::ostream& operator<<( std::ostream& out, triplet_table<DATA> &other ) {
	for( dna_alpha::alphabet_type::const_iterator i = dna_alpha::alphabet.begin(); i != dna_alpha::alphabet.end(); ++i ) {
		for( dna_alpha::alphabet_type::const_iterator j = dna_alpha::alphabet.begin(); j != dna_alpha::alphabet.end(); ++j ) {
			for( dna_alpha::alphabet_type::const_iterator k = dna_alpha::alphabet.begin(); k != dna_alpha::alphabet.end(); ++k ) {
				out << *k << *j << *i << "\t" << (double) other(*k,*j,*i) << std::endl;
			}
		}
	}
	return out;
}
};
#endif
