#ifndef BITPACKED_DNA_DATA_H_
#define BITPACKED_DNA_DATA_H_

// Haha, very hacky. Not a proper container probably. But it required fewest code changes.

#include <vector>
#include <btl/dna.h>

template< class Sequence_data, typename Alphabet >
class bitpacked_dna_data {
	private:
		typedef unsigned char byte;
		std::vector<byte> data;
	
	public:
		typedef Sequence_data sequence_data;
		typedef Alphabet alphabet;
		typedef std::vector<byte>::const_iterator const_iterator;
		
		bitpacked_dna_data() {}
		
		bitpacked_dna_data( Sequence_data const &indata ) {
			using namespace std;
			// compute the amount of space we need, we pack 4 symbols in a char
			data.resize( indata.size()/4 + 1, 0 );

			// pack the data
			unsigned int current_byte = 0;
			unsigned int symbol_num = 0;
			for( typename Sequence_data::const_iterator i = indata.begin(); i != indata.end(); i++ ) {
				byte sym = (byte) btl::dna::symbol(*i).index();
				sym <<= symbol_num;
				data[current_byte] = data[current_byte] | sym;
				
				// adjust pointers
				symbol_num += 2;
				if( symbol_num == 8 ) {
					symbol_num = 0;
					current_byte++;
				}
			}
		}
		
		// unpack n characters
		template< typename output_iterator >
		void unpack( unsigned int n, output_iterator out ) const {
			unsigned int current_byte = 0, symbol_num = 0, output_count = 0;
			while( output_count < n ) {
				byte symbol = (data[current_byte] & (3 << symbol_num)) >> symbol_num;
				*out = Alphabet::symbol( btl::dna::symbol::from_index( symbol ) );
				out++;
				output_count++;
		
				symbol_num += 2;
				if( symbol_num == 8 ) {
					symbol_num = 0;
					current_byte++;
				}
			}
		}
		
		// provide access to byte data, begin iterator
		std::vector<byte>::const_iterator begin() const {
			return data.begin();
		}
		
		// provide access to byte data, end iterator
		std::vector<byte>::const_iterator end() const {
			return data.end();
		}
		
	template < class S, typename A> friend bool operator== ( bitpacked_dna_data<S,A> const &a, bitpacked_dna_data<S,A> const &b );
	template < class S, typename A> friend std::ostream& operator<< ( std::ostream &, bitpacked_dna_data<S,A> const &b );
};

template< class Sequence_data, typename Alphabet >
inline bool operator==( bitpacked_dna_data<Sequence_data,Alphabet> const &a, bitpacked_dna_data<Sequence_data,Alphabet> const &b ) {
	return a.data == b.data;
}

template< class Sequence_data, typename Alphabet >
inline std::ostream& operator<<( std::ostream &out, bitpacked_dna_data<Sequence_data,Alphabet> const &a ) {
	typedef typename bitpacked_dna_data<Sequence_data,Alphabet>::byte byte;
	unsigned int symbol_num = 0;
	for( unsigned int current_byte = 0; current_byte < a.data.size(); ) {
		byte symbol = (a.data[current_byte] & (3 << symbol_num)) >> symbol_num;
		out << btl::dna::symbol::from_index( symbol );
		
		symbol_num += 2;
		if( symbol_num == 8 ) {
			symbol_num = 0;
			current_byte++;
		}
	}
	return out;
}

#endif /*BITPACKED_DNA_DATA_H_*/
