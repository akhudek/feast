#include "fragment_algorithm.h"

#include <iostream>
#include <vector>

namespace fragment_algorithm {

local_chain_vector_ptr positive_local_chains(	transformed_fragment_ptr_vector &fragments, double hcost, double vcost) {
	using namespace std;

	// First we sort the input fragments.
	sort(fragments.begin(), fragments.end(), tfptr_xylt());

	// Remove duplicates caused by multiple seeds.
	{
		unsigned int i, j, k;
		for( i = 1, j = 1; j < fragments.size(); ++i, ++j ) {
			if(    (*fragments[j-1])[0].a == (*fragments[j])[0].a
				&& (*fragments[j-1])[1].a == (*fragments[j])[1].a ) ++j;
			fragments[i] = fragments[j];
		}
		for( k = j-i; k > 0; --k ) fragments.pop_back();
	}

	// Now we build our index structure for easy traversal.
	vector<unsigned int> xskip;
	for (unsigned int i = 0, cx = 0; i < fragments.size(); ++i) {
		if (fragments[i]->a().x != fragments[cx]->a().x) cx = i;
		xskip.push_back(cx);
	}

	// We then do dynamic programming to find optimal chain. This is an
	// O(n^2) algorithm where n is the number of fragments.
	//
	// Each element in the array has three values. The first points to the
	// previous element in the chain, the second to the first and the second is the score.
	local_chain_vector_ptr dparray(new local_chain_vector(fragments.size()));

	int y_max = (int) (-1.0 / vcost), x_max = (int) (-1.0 / hcost);
	ereal one;
	one.set_base(1.0);

	for (unsigned int i = 0; i < fragments.size(); i++) {
		// initialize element to (-1, score of fragment i)
		ereal score = fragments[i]->score();
		(*dparray)[i].previous = -1;
		(*dparray)[i].score = score;
		(*dparray)[i].range_a = (*fragments[i])[0];
		(*dparray)[i].range_b = (*fragments[i])[1];
		(*dparray)[i].components = 1;

		// Now we walk backwards only as much as we need to. For condition handles special
		// case where there is no previous i.
		for (int j = i - 1; j >= 0; --j) {

			//cerr << i << "\t" << j << "\t" << fragments[i]->a() << "\t" << fragments[j]->a() << "\t" << x_max << "\t" << y_max << endl;

			// no need to look further
			if (fragments[i]->a().x - fragments[j]->a().x > x_max)
				break;

			// need to either binary search or jump to the correct place.
			if (abs(fragments[i]->a().y - fragments[j]->a().y) > y_max) {
				// need to binary search.
				if (fragments[i]->a().y > fragments[j]->a().y)
					j = upper_bound(fragments.begin() + xskip[j],
							fragments.begin() + j, fragments[i]->a().y - y_max,
							tfptr_ylt()) - fragments.begin();
				// otherwise we can skip right to the next i
				else
					j = xskip[j];

				// don't chain fragments that are not collinear
			} else if ((*fragments[i])[0].a <= (*fragments[j])[0].a
					|| (*fragments[i])[1].a <= (*fragments[j])[1].a) {
				continue;
				// evaluate this fragment for chaining
			} else {
				double ydiff = (double) abs(fragments[i]->a().y - fragments[j]->a().y);
				double xdiff = (double) (fragments[i]->a().x - fragments[j]->a().x);
				ereal lsh, lsv;
				lsh.set_base(hcost * xdiff);
				lsv.set_base(vcost * ydiff);
				ereal score = one * lsh * lsv;

				// if chain is better than choose it
				if (score.as_base() > 0.0 && (*dparray)[j].score * score
						> (*dparray)[i].score) {
					(*dparray)[i].previous = j;
					(*dparray)[i].score = (*dparray)[j].score * score;
					(*dparray)[i].components = (*dparray)[j].components + 1;
					(*dparray)[i].range_a.a = (*dparray)[j].range_a.a;
					(*dparray)[i].range_b.a = (*dparray)[j].range_b.a;
				}
			}
		}
	}
	return dparray;
}

bool overlap_free( local_chain_vector &chains, int start ) {
	for( int i = start; i != -1; i = chains[i].previous ) {
		if( chains[i].visited ) return false;
		chains[i].visited = true;
	}
	return true;
}


/*
 void for_each_range_no_anchor(
 closed_interval_vector &range,
 fragment_ptr_vector &chain,
 fragment_range_visitor &vis) {
 using namespace boost::lambda;

 vector<int> start, end;
 start = rangeStarts( range );
 fragment_ptr_vector_iter anchor = chain.begin();

 // adjust for no anchor regions
 if( chain.size() == 0 ) {
 end = rangeEnds( range );
 } else {
 // set to start of anchor midpoint
 end = get_anchor_midpoint( **anchor );
 }

 // do action before the first range
 vis.begin_action();

 while( 1 ) {
 // do action on anchor before following range
 if( anchor != chain.end() ) {
 // set to anchor midpoint
 end = get_anchor_midpoint( **anchor );

 vis.pre_anchor_action( **anchor );
 } else {
 end = rangeEnds( range );
 vis.end_action();
 }

 // do action on a range
 vis.range_action( start, end );

 if( anchor == chain.end() ) break;

 // do action on fragment
 vis.post_anchor_action( **anchor );

 // set next start to anchor midpoint + 1
 start = end;
 for_each( start.begin(), start.end(), _1++ );

 // move to next anchor
 anchor++;
 }
 }


 void for_each_range( closed_interval_vector &range, fragment_ptr_vector &chain, fragment_range_visitor &vis) {
 using namespace boost::lambda;

 vector<int> start, end;
 start = rangeStarts( range );
 fragment_ptr_vector_iter anchor = chain.begin();

 // adjust for no anchor regions
 if( chain.size() == 0 ) {
 end = rangeEnds( range );
 } else {
 // set to start of anchor -1
 end = (**anchor).getStart();
 for_each( end.begin(), end.end(), _1-- );
 }

 // do action before the first range
 vis.begin_action();

 while( 1 ) {
 // do action on anchor before following range
 if( anchor != chain.end() ) {
 // set to start of anchor -1
 end = (**anchor).getStart();
 for_each( end.begin(), end.end(), _1-- );

 vis.pre_anchor_action( **anchor );
 } else {
 end = rangeEnds( range );
 vis.end_action();
 }

 // do action on a range
 vis.range_action( start, end );

 if( anchor == chain.end() ) break;

 // do action on fragment
 vis.post_anchor_action( **anchor );

 // set next start to anchor end + 1
 start = (**anchor).getEnd();
 for_each( start.begin(), start.end(), _1++ );

 // move to next anchor
 anchor++;

 }
 }
 */
}
;
