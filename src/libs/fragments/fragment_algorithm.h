#ifndef __FRAGMENTALGORITHM__
#define __FRAGMENTALGORITHM__
#include <functional>
#include <limits>
#include <btl/logspace.h>

// Collection of algorithms for working with fragments.
#include "hash_hit_data.h"
#include "fast_hash_hit_data.h"
#include "fragment.h"
#include "scored_fragment.h"
#include "transformed_fragment.h"

namespace fragment_algorithm {

template<typename OutputIter>
void create_fragments(hash_hit_data &data, OutputIter out, unsigned int limit = 10000) {
	using namespace std;

	hash_hit_data::hash_data const &hits = data.get_hash();
	unsigned int numSequences = data.get_num_sequences();

	for (hash_hit_data::hash_data_citer i = hits.begin(); i != hits.end(); i++) {
		// ensure we have hits in all sequences
		bool goodhit = true;
		if (i->second.size() < numSequences) {
			goodhit = false;
		} else {
			for (unsigned int j = 0; j < numSequences; j++) {
				if (i->second[j].size() == 0)
					goodhit = false;
			}
		}

		// enumerate the hits to a set of fragments
		if (goodhit) {
			vector<vector<unsigned int> > newFragments;
			vector<vector<unsigned int> >::const_iterator j = i->second.begin();
			vector<unsigned int> numList;
			vector<unsigned int>::const_iterator k;

			// add the first sequence start and end
			for (k = j->begin(); k != j->end(); k++) {
				newFragments.push_back(vector<unsigned int> (1, *k));
			}
			j++;

			// enumerate the rest
			for (; j != i->second.end(); j++) {
				// for each fragment
				unsigned int numFrag = newFragments.size();
				for (unsigned int m = 0; m < numFrag; m++) {
					// make copies of current partial frag
					for (unsigned int n = 1; n < j->size(); n++) {
						numList = newFragments[m];
						numList.push_back((*j)[n]);
						newFragments.push_back(numList);
					}
					// for last element use existing frag
					newFragments[m].push_back((*j)[0]);
				}
			}

			// if number of fragments is less than the limit add these to the
			// result set
			// if( newFragments.size() < FRAGMENT_LIMIT ) {
			//transform( i->first.begin(), i->first.end(), ostream_iterator<char>(cout), dna_iupac_gap::decode );
			//cout << "\t" << newFragments.size() << endl;
			for (j = newFragments.begin(); j != newFragments.end(); j++) {
				scored_fragment_ptr newf(new scored_fragment());
				for (k = j->begin(); k != j->end(); k++) {
					closed_interval newinterval((int) *k, (int) *k + data.get_match_length() - 1);
					newf->add(newinterval);
				}
				// add fragment to list
				*out = newf;
				++out;
			}
				//    } else {
				//       cerr << i->first << "\t" << newFragments.size() << " filtered" << endl;
				//    }
		}
	}
}

template<typename ALPHABET, typename OutputIter>
void create_fragments(fast_hash_hit_data<ALPHABET> &data, OutputIter out ) {
	using namespace std;

	typename fast_hash_hit_data<ALPHABET>::hash_data const &hits = data.get_hash();

	vector<typename fast_hash_hit_data<ALPHABET>::position_list::const_iterator>
		position_breaks( data.number_of_sequences() + 1 ) ,
		cur_position( data.number_of_sequences() );

	// retrieve size of hit
	int match_length = (int) data.get_seed().length() - 1;

	for( typename fast_hash_hit_data<ALPHABET>::hash_data::const_iterator i = hits.begin(); i != hits.end(); i++) {
		// ensure hits in all sequences and save start positions
		bool complete = true;

		// if last sequence is empty then end matches beginning
		if( *(i->second.sequence_ends.begin()) == i->second.positions.begin() ) continue;

		// setup up pointers to sequence break points
		typename vector<typename fast_hash_hit_data<ALPHABET>::position_list::const_iterator>::iterator pb = position_breaks.begin(), cp;
		*pb = i->second.positions.begin();
		++pb;
		*pb = i->second.sequence_ends.front();
		++pb;

		// successive end points indicate empty position list
		for( typename fast_hash_hit_data<ALPHABET>::sequence_end_list::const_iterator j = i->second.sequence_ends.begin(),
			 j_next = ++(i->second.sequence_ends.begin());
			 j_next != i->second.sequence_ends.end(); ++j, ++j_next, ++pb ) {

			if( *j == *j_next ) complete = false;
			*pb = *j_next;
		}
		if( !complete ) continue;

		// complete position, now we enumerate
		for(cp = cur_position.begin(), pb = position_breaks.begin(); cp != cur_position.end(); ++cp, ++pb ) *cp = *pb;
		typename vector<typename fast_hash_hit_data<ALPHABET>::position_list::const_iterator>::reverse_iterator rpb, rcp;

		while( cur_position.front() != *(++position_breaks.begin()) ) {
			scored_fragment_ptr newf(new scored_fragment());

			for(rcp = cur_position.rbegin(); rcp != cur_position.rend(); ++rcp ) {
				closed_interval newinterval((int) **rcp, (int) **rcp + match_length);
				newf->add(newinterval);
			}

			// add fragment to list
			*out = newf;
			++out;

			// move cp pointers
			rpb = position_breaks.rbegin(), rcp = cur_position.rbegin();
			++(*rcp);
			while( *rcp == *rpb && cur_position.front() != *(++position_breaks.begin()) ) {
				++rpb;
				*rcp = *rpb;
				++rcp;
				++(*rcp);
			}
		}
	}
}

struct greater_fragment: public std::binary_function<scored_fragment_ptr,
		scored_fragment_ptr, bool> {
	inline bool operator()(scored_fragment_ptr a, scored_fragment_ptr b) {
		return *a >> *b;
	}
};

struct null_cf: public std::binary_function<scored_fragment, scored_fragment,
		double> {
	inline ereal operator()(scored_fragment const &J, scored_fragment const &I) {
		return ereal(1.0);
	}
};

struct cost_bonus_cf: public std::binary_function<scored_fragment,
		scored_fragment, double> {
	inline ereal operator()(scored_fragment const &J, scored_fragment const &I) {
		// J >> I
		int x1 = J[1].b, y1 = J[0].b, x2 = I[1].a, y2 = I[0].a;
		ereal r;

		r.set_base((double) abs(y2 - y1) - (double) abs(x1 - x2 + y2 - y1));
		return r;
	}
};

struct gap_score_cf: public std::binary_function<scored_fragment,
		scored_fragment, double> {
	ereal cost;
	gap_score_cf(double c) :
		cost(c) {
	}
	inline ereal operator()(scored_fragment const &J, scored_fragment const &I) {
		// J >> I
		int x1 = J[1].b, y1 = J[0].b, x2 = I[1].a, y2 = I[0].a;
		//return (double)abs(y2-y1) - (double)abs(x1-x2+y2-y1);
		ereal r;
		r.set_base(cost.as_base() * (double) abs(x1 - x2 + y2 - y1));
		return r;
	}
};

struct gap_missing_cf: public std::binary_function<scored_fragment,
		scored_fragment, double> {

	ereal gap_cost;
	ereal missing_cost;

	gap_missing_cf(double c1, double c2) :
		gap_cost(c1), missing_cost(c2) {
	}

	inline ereal operator()(scored_fragment const &J, scored_fragment const &I) {
		// J >> I
		int x1 = J[1].b, y1 = J[0].b, x2 = I[1].a, y2 = I[0].a;
		//return (double)abs(y2-y1) - (double)abs(x1-x2+y2-y1);
		ereal r;
		r.set_base(gap_cost.as_base() * (double) abs(x1 - x2 + y2 - y1)
				+ missing_cost.as_base() * (double) (y2 - y1));
		return r;
	}
};

// This algorithm computes the optimal fragment chain using a heuristic
// that is no worse than O(n^2) and hopefully much faster in most cases.
template<class OutputIterator>
void optimal_chain_fast(scored_fragment_ptr_vector fragments,
		OutputIterator output) {
	using namespace std;

	// This map stores the scores of existing chains in increasing order.
	// Each score contains the postiion of the last fragment in the
	// chain.
	typedef multimap<ereal, int> score_map;
	typedef score_map::reverse_iterator score_map_rev_iter;
	score_map current_chains;

	// First we must sort the input fragments. We will use lexicographic
	// order for this.
	sort(fragments.begin(), fragments.end(), greater_fragment());

	// We then do dynamic programming to find optimal chain. This is an
	// O(n^2) algorithm where n is the number of fragments.
	//
	// Each element in the array has two values. The first points to the
	// previous element in the chain and the second is the score.
	vector<pair<int, ereal > > dparray(fragments.size());

	// track best chain
	pair<int, ereal > best_chain;
	best_chain.first = -1;
	best_chain.second = 0.0; // we can have negative good anchors?

	for (unsigned int i = 0; i < fragments.size(); i++) {
		// initialize element to (-1, score of fragment i)
		ereal score = fragments[i]->score();
		dparray[i].first = -1;
		dparray[i].second = score;

		// Find max scoring previous chain that does not conflict with
		// this fragment.
		for (score_map_rev_iter c = current_chains.rbegin(); c
				!= current_chains.rend(); c++) {
			// if the fragment does not conflict, then we add it and quite
			// because it is the highest non-conflicting fragment.

			// If j >> i for fagments i and j then we can chain.
			if (*(fragments[c->second]) >> *(fragments[i])) {
				// chain is better, choose it and exit for
				dparray[i].first = c->second;
				dparray[i].second = c->first * score;
				break; // short circuit rest of search
			}
		}

		// add current fragment to map
		current_chains.insert(make_pair(dparray[i].second, i));

		// check for best chain
		if (dparray[i].second > best_chain.second) {
			best_chain.first = i;
			best_chain.second = dparray[i].second;
		}
	}

	for (int i = best_chain.first; i != -1; i = dparray[i].first, ++output) {
		*output = fragments[i];
	}
}

// This algorithm computes the optimal fragment chain
template<typename ChainCostFn, typename OutputIterator>
void optimal_chain(scored_fragment_ptr_vector fragments, OutputIterator output,
		ChainCostFn chain_score_fn) {
	using namespace std;

	// First we must sort the input fragments. We will use lexicographic
	// order for this.
	sort(fragments.begin(), fragments.end(), greater_fragment());

	// We then do dynamic programming to find optimal chain. This is an
	// O(n^2) algorithm where n is the number of fragments.
	//
	// Each element in the array has two values. The first points to the
	// previous element in the chain and the second is the score.
	vector<pair<int, ereal > > dparray(fragments.size());

	// track best chain
	pair<int, ereal > best_chain;
	best_chain.first = -1;
	best_chain.second = 0.0;

	for (unsigned int i = 0; i < fragments.size(); i++) {
		// initialize element to (-1, score of fragment i)
		ereal score = fragments[i]->score();
		dparray[i].first = -1;
		dparray[i].second = score;

		for (unsigned int j = 0; j < i; j++) {
			// If j >> i for fagments i and j then we can chain.
			if (*fragments[j] >> *fragments[i]) {
				//cerr << *fragments[j] << "\t" << *fragments[i] << endl;
				ereal chain_score = chain_score_fn( *fragments[j], *fragments[i]);
				// if chain is better than choose it
				if (dparray[j].second * chain_score * score > dparray[i].second) {
					dparray[i].first = j;
					dparray[i].second = dparray[j].second * chain_score * score;
				}
			}
		}

		// check for best chain
		if (dparray[i].second > best_chain.second) {
			best_chain.first = i;
			best_chain.second = dparray[i].second;
		}
	}
	for (int i = best_chain.first; i != -1; i = dparray[i].first, ++output) {
		*output = fragments[i];
	}
}

struct local_chain {
	ereal score;
	closed_interval range_a;
	closed_interval range_b;
	int components;
	int previous;
	bool visited;
	local_chain() : score(0.0), range_a(-1, -1), range_b(-1, -1), components(0), previous(-1), visited(false) {}
};

typedef std::vector<local_chain> local_chain_vector;
typedef boost::shared_ptr<local_chain_vector> local_chain_vector_ptr;


struct local_chain_ltscore {
	bool operator()(local_chain const &a, local_chain const &b) {
		return a.score < b.score;
	}
};

struct local_chain_id_gtscore {
	local_chain_vector &chain;
	local_chain_id_gtscore( local_chain_vector &chains ) : chain(chains) {}
	bool operator()(int a, int b) { return chain[a].score > chain[b].score; }
};

struct local_chain_eqstart {
	bool operator()(local_chain const &a, local_chain const &b) {
		return a.range_a.a == b.range_a.a && a.range_b.a == b.range_b.a;
	}
};


// This algorithm computes the optimal fragment chain
template<typename ChainCostFn>
local_chain_vector_ptr local_chains(scored_fragment_ptr_vector fragments, ChainCostFn chain_score_fn) {
	using namespace std;

	// First we must sort the input fragments. We will use lexicographic
	// order for this.
	sort(fragments.begin(), fragments.end(), greater_fragment());

	// We then do dynamic programming to find optimal chain. This is an
	// O(n^2) algorithm where n is the number of fragments.
	//
	// Each element in the array has three values. The first points to the
	// previous element in the chain, the second to the first and the second is the score.
	local_chain_vector_ptr dparray(new local_chain_vector(fragments.size()));

	for (unsigned int i = 0; i < fragments.size(); i++) {
		// initialize element to (-1, score of fragment i)
		ereal score = fragments[i]->score();
		(*dparray)[i].previous = -1;
		(*dparray)[i].score = score;
		(*dparray)[i].range_a = (*fragments[i])[0];
		(*dparray)[i].range_b = (*fragments[i])[1];
		(*dparray)[i].components = 1;

		for (unsigned int j = 0; j < i; j++) {
			// If j >> i for fagments i and j then we can chain.
			if (*fragments[j] >> *fragments[i]) {
				//cerr << *fragments[j] << "\t" << *fragments[i] << endl;
				ereal chain_score = chain_score_fn(
						*fragments[j], *fragments[i]);
				// if chain is better than choose it
				if ((*dparray)[j].score * chain_score * score
						> (*dparray)[i].score) {
					(*dparray)[i].previous = j;
					(*dparray)[i].score = (*dparray)[j].score * chain_score
							* score;
					(*dparray)[i].components = (*dparray)[j].components + 1;
					(*dparray)[i].range_a.a = (*dparray)[j].range_a.a;
					(*dparray)[i].range_b.a = (*dparray)[j].range_b.a;
				}
			}
		}
	}
	sort(dparray->begin(), dparray->end(), local_chain_ltscore());
	return dparray;
}

bool overlap_free( local_chain_vector &chains, int start );


struct tfptr_yxlt: public std::binary_function<transformed_fragment_ptr,
		transformed_fragment_ptr, bool> {
	tfptr_yxlt() {
	}

	inline bool operator()(transformed_fragment_ptr const &I,
			transformed_fragment_ptr const &J) {
		if (I->a().y < J->a().y)
			return true;
		if (I->a().y == J->a().y && I->a().x < J->a().x)
			return true;
		return false;
	}
};

struct tfptr_xylt: public std::binary_function<transformed_fragment_ptr,
		transformed_fragment_ptr, bool> {
	tfptr_xylt() {
	}

	inline bool operator()(transformed_fragment_ptr const &I,
			transformed_fragment_ptr const &J) {
		if (I->a().x < J->a().x)
			return true;
		if (I->a().x == J->a().x && I->a().y < J->a().y)
			return true;
		return false;
	}
};

struct tfptr_ylt: public std::binary_function<int const&,
		transformed_fragment_ptr const&, bool> {
	tfptr_ylt() {
	}
	inline bool operator()(int const &v, transformed_fragment_ptr const &I) {
		return I->a().y < v;
	}
};

// Precondition is that fragments all have score 1.0.
local_chain_vector_ptr positive_local_chains(	transformed_fragment_ptr_vector &fragments, double hcost, double vcost);

inline std::pair<int, int> midpoint(fragment const &f) {
	assert( f.size() == 2 );
	assert( f[0].size() == f[1].size() );
	int offset = (int) f[0].size() / 2;
	return std::make_pair(f[0].a + offset, f[1].a + offset);
}

/*
 fragment_ptr_vector_ptr optimal_chain( fragment_ptr_vector &frags );

 void for_each_range_no_anchor( closed_interval_vector &range, fragment_ptr_vector &chain, fragment_range_visitor &vis );

 void for_each_range( closed_interval_vector &range, fragment_ptr_vector &chain, fragment_range_visitor &vis );

 // assign a score to each fragment
 template<typename ScoreFunction> inline void assign_scores( fragment_ptr_vector &frags, ScoreFunction &scoreFunc ) {
 for( fragment_ptr_vector_iter i = frags.begin(); i != frags.end(); i++ ) {
 (**i).setScore( scoreFunc( **i ) );
 //cerr << **i << endl;
 }
 }

 */

}
;
#endif
