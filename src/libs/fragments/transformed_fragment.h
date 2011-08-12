#ifndef TRANSFORMED_FRAGMENT_H__
#define TRANSFORMED_FRAGMENT_H__

#include "ereal.h"

#include "fragment.h"
#include "point.h"
#include <btl/logspace.h>
#include <boost/shared_ptr.hpp>

// Only supports 2d fragments.
class transformed_fragment : public fragment {

private:
	double duplicity;
	ereal sc;
	point pa;
	point pb;

	static double const PI4;
public:
	transformed_fragment( fragment const &f );

	// set the score
	void set_score( ereal s );

	// get the score
	ereal score() const;

	// get transformed coordinate
	point const &a() const;
	point const &b() const;
};

inline transformed_fragment::transformed_fragment( fragment const &f ) : fragment(f) {
	assert( f.size() == 2 );
	using namespace std;
	pa = point(f[0].a,f[1].a);
	pa.rotate(-PI4);
	pb = point(f[0].b,f[1].b);
	pb.rotate(-PI4);
}

inline point const & transformed_fragment::a() const { return pa; }
inline point const & transformed_fragment::b() const { return pb; }


inline void transformed_fragment::set_score( ereal s ) { sc = s; }

inline ereal transformed_fragment::score() const { return sc; }

// define some shared pointers
typedef boost::shared_ptr<transformed_fragment> transformed_fragment_ptr;

// define fragment vectors
typedef std::vector<transformed_fragment_ptr> transformed_fragment_ptr_vector;
typedef transformed_fragment_ptr_vector::iterator transformed_fragment_ptr_vector_iter;
typedef transformed_fragment_ptr_vector::reverse_iterator transformed_fragment_ptr_vector_riter;
typedef transformed_fragment_ptr_vector::const_iterator transformed_fragment_ptr_vector_citer;
typedef boost::shared_ptr<transformed_fragment_ptr_vector> transformed_fragment_ptr_vector_ptr;

// define transformed_fragment lists
typedef std::list<transformed_fragment_ptr> transformed_fragment_ptr_list;
typedef transformed_fragment_ptr_list::iterator transformed_fragment_ptr_list_iter;
typedef transformed_fragment_ptr_list::const_iterator transformed_fragment_ptr_list_citer;
typedef boost::shared_ptr<transformed_fragment_ptr_list> transformed_fragment_ptr_list_ptr;

#endif

