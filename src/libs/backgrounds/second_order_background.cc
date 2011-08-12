#include "second_order_background.h"

second_order_background::second_order_background() {}

// convert each value to log base sqrt(2)
void second_order_background::set_conditionals( triplet_table<logint> &newq ) {
   cond = newq;
}

second_order_background::preprocess_data_ptr second_order_background::preprocess( dna_sequence_data &a, closed_interval r ) {
   // some obvious assertions
   assert( r.a <= r.b );
   assert( r.a >= 0 ); 
   assert( r.b < (int) a.size() );
   assert( r.size() >= 3 );

   // create data structure, reserve memory, and save the range
   preprocess_data_ptr ppd( new preprocess_data() );
   ppd->data.reserve( r.size() );
   ppd->range = r;
   
   // First two positions we set to 1 since we can't do much else. I suppose we base probabilities over seeing pairs, but this 
   // just makes things too complicated and likely isn't worth it.
   logint last = 1.0;
   ppd->data.push_back( last ); ppd->data.push_back( last );
   
   // start the computation
   for( int i = r.a+2; i <= r.b; i++ ) {
      last *= cond(a[i],a[i-1],a[i-2]);
      ppd->data.push_back(last);
   }
   return ppd;
}



