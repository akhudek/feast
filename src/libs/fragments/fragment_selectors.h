#ifndef FRAGMENT_SELECTORS_H__
#define FRAGMENT_SELECTORS_H__
#include <functional>
//#include "fpint.h"

/*struct fragment_score : public std::unary_function<scored_fragment_ptr,fpint> {
   public:
	inline ereal operator() ( scored_fragment_ptr f ) { return f->score(); }
};*/

struct by_increasing_fragment_score : public std::binary_function<scored_fragment_ptr,scored_fragment_ptr,bool> {
   public:
      inline bool operator() ( scored_fragment_ptr a, scored_fragment_ptr b ) {
         return a->score() < b->score();
      }
};

/*
struct by_increasing_extension_score : public std::binary_function<scored_fragment_ptr,scored_fragment_ptr,bool> {
public:
	inline bool operator() ( scored_fragment_ptr a, scored_fragment_ptr b ) {
		int aext = a->forward_extension() + a->backward_extension();
		int bext = b->forward_extension() + b->backward_extension();
		if( aext == bext ) return a->score() < b->score();
		return aext < bext;
	}
};
*/



#endif
