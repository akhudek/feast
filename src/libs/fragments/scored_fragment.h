#ifndef SCORED_FRAGMENT_H__
#define SCORED_FRAGMENT_H__

#include "fragment.h"
#include "local_segmentation.h"


class scored_fragment : public fragment {
public:
	typedef local_segmentation<>::result score_data;
	unsigned int pset;

private:
	score_data pscore;

public:
	scored_fragment() {}

	// set the score
	void set_score( ereal s );

	// get the score
	ereal score() const;

	// set the score data
	void set_score_data( score_data const &sd);

	// set region i type
	void set_parameter_set( unsigned int s );

	// get region i
	unsigned int parameter_set( ) const;
};

inline void scored_fragment::set_score( ereal s ) { pscore.score = s; }

inline void scored_fragment::set_score_data( score_data const &sd ) { pscore = sd; }

inline ereal scored_fragment::score() const { return pscore.score; }

inline void scored_fragment::set_parameter_set( unsigned int p ) { pset = p; }

inline unsigned int scored_fragment::parameter_set() const { return pset; }

// define some shared pointers
typedef boost::shared_ptr<scored_fragment> scored_fragment_ptr;

// define fragment vectors
typedef std::vector<scored_fragment_ptr> scored_fragment_ptr_vector;
typedef scored_fragment_ptr_vector::iterator scored_fragment_ptr_vector_iter;
typedef scored_fragment_ptr_vector::reverse_iterator scored_fragment_ptr_vector_riter;
typedef scored_fragment_ptr_vector::const_iterator scored_fragment_ptr_vector_citer;
typedef boost::shared_ptr<scored_fragment_ptr_vector> scored_fragment_ptr_vector_ptr;

// define scored_fragment lists
typedef std::list<scored_fragment_ptr> scored_fragment_ptr_list;
typedef scored_fragment_ptr_list::iterator scored_fragment_ptr_list_iter;
typedef scored_fragment_ptr_list::const_iterator scored_fragment_ptr_list_citer;
typedef boost::shared_ptr<scored_fragment_ptr_list> scored_fragment_ptr_list_ptr;

#endif

