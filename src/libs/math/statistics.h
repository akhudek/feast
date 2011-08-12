#ifndef STATISTICS_H__
#define STATISTICS_H__
#include <cmath>

namespace statistics {
   // take the mean
   template< typename Iter, typename Selector > 
   typename Selector::result_type mean( Iter start, Iter end, Selector value ) {
      typename Selector::result_type sum = 0, n = 0;
      for( Iter i = start; i != end; i++, n++ ) sum += value(*i);
      return sum/n;
   }

   // Sample standard deviation, this just converts everthing to double do to the computation. 
   // Probably could be much better, but I need it to work sooner than later.
   template< typename Iter, typename Selector > 
   typename Selector::result_type std_deviation( Iter start, Iter end, Selector value ) {
      double m = (double)mean( start, end, value ), n = 0, sum = 0, v;
      for( Iter i = start; i != end; i++, n++ ) {
         v = (double)value(*i) - m;
         sum += v*v;
      }
      typename Selector::result_type result( std::sqrt( sum/(n-1.0) ));
      return result;
   }

};
#endif
