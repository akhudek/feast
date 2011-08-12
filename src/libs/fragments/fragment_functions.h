#ifndef FRAGMENT_FUNCTIONS_H__
#define FRAGMENT_FUNCTIONS_H__
#include <functional>
#include "dna_sequence.h"

// given a fragment, build two fragments, one for each end
template<typename FRAGMENT>
inline std::pair<boost::shared_ptr<FRAGMENT>,boost::shared_ptr<FRAGMENT> > build_ends( FRAGMENT &f ) {
	boost::shared_ptr<FRAGMENT> left( new FRAGMENT() ), right( new FRAGMENT() );
	left->add( closed_interval(f[0].a-1,f[0].a-1) );
	left->add( closed_interval(f[1].a-1,f[1].a-1) );
	right->add( closed_interval(f[0].b+1,f[0].b+1) );
	right->add( closed_interval(f[1].b+1,f[1].b+1) );
	return std::make_pair(left,right);
}

// Return true if the score is below threshold t.
struct score_below : public std::unary_function<scored_fragment_ptr,bool> {
   private:
	ereal t;
   public:
	score_below( ereal new_t ) : t(new_t) {}
      inline bool operator()( scored_fragment_ptr f ) const { return f->score() < t; }
};

/*
// Return true if either forward or backward extension are below t
struct extension_below : public std::unary_function<scored_fragment_ptr,bool> {
   private:
      int e;
   public:
      extension_below( int new_e ) : e(new_e) {}
      inline bool operator()( scored_fragment_ptr f ) const {return (f->forward_extension() < e) || (f->backward_extension() < e); }
};*/


// Score fragment.
template<typename LocalAligner>
struct assign_score : public std::unary_function<scored_fragment_ptr,void> {
   private:
      LocalAligner &aligner;
      dna_sequence &seqa;
      dna_sequence &seqb;
      int min_ext;

   public:
      assign_score( LocalAligner &aa, dna_sequence &a, dna_sequence &b, int mext = 0 )
         : aligner(aa) , seqa(a), seqb(b), min_ext(mext) {}
      inline void operator()( scored_fragment_ptr f ) {
         using namespace boost;
         typename LocalAligner::result rf = aligner.extend_forward( seqa, seqb, (*f)[0].a, (*f)[1].a, min_ext );
         //typename LocalAligner::result rb = aligner.extend_backward( seqa, seqb,(*f)[0].b, (*f)[1].b, min_ext );
         (*f).set_score( rf.score );
         //(*f).set_forward( rf );
         //(*f).set_backward( rb );
      }
};

// Return a fragment representing the region between two fragments. End points are used as boundaries.
template< typename FRAGMENT >
inline FRAGMENT interior( FRAGMENT const &x, FRAGMENT const &y ) {
   FRAGMENT new_fragment;
   new_fragment.add(closed_interval(x[0].b+1, y[0].a-1));
   new_fragment.add(closed_interval(x[1].b+1, y[1].a-1));
   return new_fragment;
}


// Compute the fragment representing the region between the two fragments. Mid point and Mid point - 1 are
// boundaries.
template<typename FRAGMENT>
inline FRAGMENT midpoint_interior( FRAGMENT const &x, FRAGMENT const &y ) {
   int mid_x = (int)x[0].size()/2;
   int mid_y = (int)y[0].size()/2;
   FRAGMENT new_fragment;
   new_fragment.add(closed_interval(x[0].a+mid_x,y[0].a+mid_y-1));
   new_fragment.add(closed_interval(x[1].a+mid_x,y[1].a+mid_y-1));
   return new_fragment;
}

// Compute the area of the region between two fragments.
inline double interior_area( fragment const &x, fragment const &y ) {
   if( y[0].a == x[0].b+1 ) return 0;
   if( y[1].a == x[1].b+1 ) return 0;
   int first_size = (y[0].a-1) - (x[0].b+1) + 1;
   int second_size = (y[1].a-1) - (x[1].b+1) + 1;
   return (double)first_size*(double)second_size;
}

template <typename SEQDATA>
inline void expand_closed_interval( SEQDATA &data, closed_interval &c, int min_size ) {
   if( c.size() < min_size ) {
      // expand each side by an equal amount
      int expand_right = 0, expand_left = 0;
      unsigned int missing = min_size - c.size();
      expand_right = expand_left = missing/2;

      if( c.a - expand_left < 0 ) {
         // we collide left
         expand_left = c.a;
         expand_right = missing - expand_left;

      } else if( c.b + expand_right >= (int)data.size() ) {
         // we collide right
         expand_right = data.size() - c.b - 1;
         expand_left = missing - expand_right;
      }
      c.a -= expand_left;
      c.b += expand_right;
   }
	if( c.a < 0 ) c.a = 0;
   if( c.b >= (int)data.size() ) c.b = (int)data.size()-1;
}



#endif

