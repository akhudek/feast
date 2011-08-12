/** Positive local chaining algorithm.
 */

#ifndef CHAINING_ALGORITHM_H_
#define CHAINING_ALGORITHM_H_

#include "ereal.h"
#include <boost/shared_ptr.hpp>
#include <vector>
#include <algorithm>
#include <functional>
#include "local_chain.h"
#include "algorithm_extra.h"
#include "point.h"

template< typename PointIterator >
boost::shared_ptr< std::vector<local_chain> > positive_local_chains( PointIterator pit_start, PointIterator pit_end, double hcost, double vcost) {
	using namespace std;
	using namespace boost;
	using namespace btl;

	// Transform points.
	shared_ptr< vector<local_chain> > dparray( new vector<local_chain>() );
	for( PointIterator i = pit_start; i != pit_end; ++i ) dparray->push_back( local_chain(*i) );

	// First we sort the input fragments.
	sort(dparray->begin(), dparray->end(), local_chain_xylt());

	// Remove duplicates caused by multiple seeds.
	{
		unsigned int i, j, k;
		for( i = 1, j = 1; j < dparray->size(); ++i, ++j ) {
			if( (*dparray)[j-1].p == (*dparray)[j].p ) ++j;
			(*dparray)[i] = (*dparray)[j];
		}
		for( k = j-i; k > 0; --k ) dparray->pop_back();
	}

	// Now we fill in our index structure for easy traversal.
		for( vector<local_chain>::iterator i = dparray->begin(), cx = dparray->begin(); i != dparray->end(); ++i ) {
		if( i->tp.x != cx->tp.x ) {
			cx = i;
			--cx;
		}
		i->xskip = cx;
	}

	// We then do dynamic programming to find optimal chain. This is an
	// O(n^2) algorithm where n is the number of fragments.
	//
	// Each element in the array has three values. The first points to the
	// previous element in the chain, the second to the first and the second is the score.

	int y_max = (int) (-1.0 / vcost), x_max = (int) (-1.0 / hcost), xdiff, ydiff;
	ereal one; one.set_base(1.0);

	dparray->front().previous = dparray->end();
	dparray->front().score = one;
	for( vector<local_chain>::iterator i = ++(dparray->begin()); i != dparray->end(); ++i ) {
		// initialize element to (-1, score of fragment i)
		i->previous = dparray->end();
		i->score = one;

		// Now we walk backwards only as much as we need to. For condition handles special
		// case where there is no previous i.
		for( vector<local_chain>::iterator j = i - 1; j != dparray->begin() - 1; --j ) {

			// no need to look further
			xdiff = i->tp.x - j->tp.x;
			if( xdiff > x_max) break;

			ydiff = i->tp.y - j->tp.y;

			// need to either binary search or jump to the correct place.
			if (abs(ydiff) > y_max) {
				// need to binary search.
				if( i->tp.y > j->tp.y )
					j = upper_bound(j->xskip, j, i->tp.y - y_max, local_chain_ylt() );
				// otherwise we can skip right to the next i
				else
					j = j->xskip;

				// don't chain fragments that are not collinear
			} else if( i->p.x <= j->p.x	|| i->p.y <= j->p.y) {
				continue;
				// evaluate this fragment for chaining
			} else {
				ereal lsh, lsv;
				lsh.set_base(hcost * (double)xdiff);
				lsv.set_base(vcost * (double)abs(ydiff) );
				ereal score = j->score * one * lsh * lsv;

				// if chain is better than choose it
				if( score > i->score ) {
					i->previous = j;
					i->score = score;
				}
			}
		}
	}
	return dparray;
}

inline bool paint_chain( std::vector<local_chain>::iterator j, std::vector<local_chain>::iterator end ) {
	using namespace std;
	for( ; j != end; j = j->previous ) {
		if( j->seen ) return false;
		j->seen = true;
	}
	return true;
}

/** Prune and filter the chain tree.
 *  We output the positions of the filtered chains to oit.
 */
template< typename OutputIterator >
void prune_and_filter_chains( std::vector<local_chain> &chains, ereal T, OutputIterator oit ) {
	using namespace std;
	using namespace btl;
	using namespace __gnu_cxx;

	// First we make a list of position score pairs with scores above T.
	typedef ereal scoret;
	typedef pair<vector<local_chain>::iterator,scoret>  ps_pair;
	vector<ps_pair> filtered_positions;

	for( vector<local_chain>::iterator i = chains.begin(); i != chains.end(); ++i ) {
		if( i->score >= T ) filtered_positions.push_back( make_pair( i, i->score ) );
	}

	// Sort filtered positions by score.
	sort( filtered_positions.begin(), filtered_positions.end(), gte_by2nd<ps_pair>() );

	// Paint and filter chains.
	for( vector<ps_pair>::const_iterator i = filtered_positions.begin(); i != filtered_positions.end(); ++i ) {
		if( paint_chain(i->first, chains.end() ) ) {
			*oit = i->first;
			++oit;
		}
	}
}

point get_chain_middle( std::vector<local_chain>::iterator i, std::vector<local_chain>::iterator end );

std::vector<local_chain>::iterator get_chain_middle_it( std::vector<local_chain>::iterator i, std::vector<local_chain>::iterator end );


template< typename OutputIterator >
void get_chain_points( std::vector<local_chain>::iterator i, std::vector<local_chain>::iterator end, OutputIterator oit ) {
	using namespace std;
	while( i != end ) {
		*oit = i->p;
		++oit;
		// Move forward in i.
		i = i-> previous;
	}
}

template< typename OutputIterator >
void get_chain_points_prioritize_middle( std::vector<local_chain>::iterator begin, std::vector<local_chain>::iterator end, OutputIterator oit ) {
	using namespace std;

	vector<local_chain>::iterator mid = get_chain_middle_it(begin,end);
	vector<point> right_half;
	vector<local_chain>::iterator i = begin;

	while( i != mid ) {
		right_half.push_back( i->p );
		i = i->previous;
	}

	vector<point>::reverse_iterator j = right_half.rbegin();
	i = mid;
	while( i != end || j != right_half.rend() ) {
		if( i != end ) {
			*oit = i->p;
			++oit;
			i = i->previous;
		}
		if( j != right_half.rend() ) {
			*oit = *j;
			++oit;
			++j;
		}
	}
}


#endif /* CHAINING_ALGORITHM_H_ */
