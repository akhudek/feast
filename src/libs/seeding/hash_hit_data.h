#ifndef HASH_HIT_DATA_H__
#define HASH_HIT_DATA_H__
// This class and hash_count_data do not use inheritance due to performance. Virtual
// functions require an extra call.
#include <vector>
#include <hash_map>
#include "dna_sequence.h"
#include "char_hash_function.h"
#include "bitpacked_dna_data.h"

class hash_hit_data {
public:
	typedef bitpacked_dna_data<dna_sequence_data, dna_sequence::alphabet>
			packed_data;
	typedef __gnu_cxx ::hash_map<packed_data, std::vector<std::vector< unsigned int> >, hashf<packed_data> > hash_data;
	typedef hash_data::iterator hash_data_iter;
	typedef hash_data::const_iterator hash_data_citer;

protected:
	hash_data hash;
	unsigned int match_length;
	unsigned int num_seqs;

public:
	hash_hit_data();

	hash_hit_data(unsigned int ml, unsigned int ns);

	unsigned int get_num_sequences() const;

	unsigned int get_match_length() const;

	inline hash_data const &get_hash() const {
		return hash;
	}

	void add(unsigned int seqi, unsigned int seqp,
			dna_sequence_data const &pattern);
};

//------------------------------------
// Inline functions.

inline unsigned int hash_hit_data::get_num_sequences() const {
	return num_seqs;
}

inline unsigned int hash_hit_data::get_match_length() const {
	return match_length;
}

inline void hash_hit_data::add(unsigned int seqi, unsigned int seq_pos,
		dna_sequence_data const &pattern) {
	packed_data packed(pattern);
	std::vector<std::vector<unsigned int> > &hash_bin = hash[packed];
	hash_bin.resize(seqi + 1);
	hash_bin[seqi].push_back(seq_pos);
}

//------------------------------------
// Some typedefs.

// pointer definition
typedef boost::shared_ptr<hash_hit_data> hash_hit_data_ptr;

// simple allocator
hash_hit_data_ptr new_hash_hit_data(unsigned int i, unsigned int j);

#endif /*HASH_HIT_DATA_H_*/
