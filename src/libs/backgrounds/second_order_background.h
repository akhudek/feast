#ifndef SECOND_ORDER_BACKGROUND_H__
#define SECOND_ORDER_BACKGROUND_H__
#include <vector>
#include <boost/shared_ptr.hpp>
#include <blitz/tinyvec.h>
#include "logint.h"
#include "closed_interval.h"
#include "dna_sequence.h"
#include "triplet_table.h"

class second_order_background {
   protected:
      // conditional probabilities
      triplet_table<logint> cond;

      // structure storing preprocess data
      struct preprocess_data {
         std::vector<logint> data;
         closed_interval range;
      };

   public:
      typedef boost::shared_ptr<preprocess_data> preprocess_data_ptr;

      second_order_background();   

      // set background 
      void set_conditionals( triplet_table<logint> &newq );

      // preprocess sequence data for fast queries
      preprocess_data_ptr preprocess( dna_sequence_data &a, closed_interval c );

      // compute model probability for a subinterval given preprocessing data
      inline logint pr( preprocess_data &ppd, int a, int b );
};

inline logint second_order_background::pr( second_order_background::preprocess_data &ppd, int a, int b ) {
   //assert( r.a >= ppd.range.a && r.a <= ppd.range.b );
   //assert( r.b >= ppd.range.a && r.b <= ppd.range.b );
   //assert( r.a <= r.b );

   logint prefix = ( a > ppd.range.a ) ? ppd.data[a - 1 - ppd.range.a] : logint(1.0);
   return ppd.data[b - ppd.range.a]/prefix;
}

#endif

