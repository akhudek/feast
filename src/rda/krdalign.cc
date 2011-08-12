#include "krdalign.h"

#include <cmath>
#include <cassert>
#include <algorithm>

//#define DEBUG_TRACEBACK
//#define DEBUG_DP 2

// initialize storage
void krdalign::init_storage() {
   Sd.resize(k+shape_z+2,k+1);
   Sr.resize(k+shape_z+2,k+1);
//   St.resize(k+1,k+1); // this is for debug only
   I.resize(k+shape_z+2,k+1);
   M.resize(k+shape_z+2,k+1);
   odds_a.resize(k+shape_z+2);
   odds_b.resize(k+shape_z+2);
}

void krdalign::init_shape() {
   using namespace std;
   double dk = (double) k;
   double dy = (double) shape_y;
   shape_z = (int)((dy + 1)*dk*dk/(dk+2.0*dy));

   // precompute the whole boundary, this saves us a fair amount of time in the long run
   // since we compute this at nearly every matrix location
   shape_fx.resize(shape_z+1);
   shape_fy.resize(k+2); // inverse of shape_fx
  
   // init accounting for shape_fy
   int last_fx = k;
   for( int i = 0; i <= (int)k; i++ ) {
      shape_fy(i) = 0;
   }
   shape_fy(last_fx) = 0;

   for( int i = 0; i <= shape_z; i++ ) {
      shape_fx(i) = shape_f(i);
   
      // deal with shape_fy   
      if( shape_fx(i) != last_fx ) {
         shape_fy(shape_fx(i)) = shape_fy(last_fx);
         last_fx = shape_fx(i);
         shape_fy(last_fx) += 1;
      } else {
         shape_fy(last_fx) += 1;
      }
   }
  
   // fill in missing values 
   for( int i = k; i >= 1; i-- ) {
      if( shape_fy(i) == 0 ) shape_fy(i) = shape_fy(i+1);
   }
}


// We assume the BLASTN default scoring scheme with Juke's/Cantor background. All logs are in base sqrt(2).
// Defaults          Adjusted to remove/include the odds portion.
// Match:       1    PrMatch:       log(4^{1}/16)    = -4     remove odds     
// Mismatch:   -3    PrMismatch:    log(4^{-3}/16)   = -20    remove odds    
// Gap extend: -2    PrGapExtend:   log(4^{-2})      = -8     odds cancel, no need to remove
//                   PrNoGapExtend: log(1-4^{-2})    = -0.19
// Gap open:   -5  combination of gap open and gap close
// 4^{-5}      = pr_gap_open*pr_gap_close 
//             = pr_gap_open*(1-pr_gap_extend)
// pr_gap_open = 4^{-5}/(1-4^{-2}) 
//                   PrGapOpen:     log(pr_gap_open)  = -19.81
//                   PrNoGapOpen:   log(1-pr_gap_open)= -0.00  too small for our precision
// Pr A,T,G,C: 1/4   Pr A,T,G,C:    log(1/4)          = -4     add up odds separately
void krdalign::init_pr() {
   // match/mismatch probabilities
   ereal ma = 0.25; // match
   ereal mm = std::pow( std::sqrt(2.0), -20.0 );
   p =   ma, mm, mm, mm,
         mm, ma, mm, mm,
         mm, mm, ma, mm,
         mm, mm, mm, ma;

   // background probabilities
   q =   ma, ma, ma, ma;

   // precompute transition probabilities
   set_popen( 0.001043033 );
   set_mean_gap_length(1.07);
}


// Set all parameters.
void krdalign::set_parameters( nw_model_parameters &mp ) {
	set_popen( mp.pr_open );
	set_mean_gap_length( mp.mean_gap_length );
	set_p( mp.p );
	set_q( mp.q );
}

// set mean gap length to n
void krdalign::set_mean_gap_length( double n ) {
   // save value in double precision
   original_b = 1.0 - 1.0/n;
 
   // Since a gap is length 1 by the virtue of opening it, we adjust n down by 1. 
   // The remaining gap extension probability corresponds to the mean of a 
   // geometric distribution allowing zero extensions.
   p_b = original_b;
   p_1mb = 1.0 - original_b;

  compute_sd();
}

// convert to log base sqrt(2)
void krdalign::set_popen( double a ) {
   using namespace std;
   p_a = a;
   p_1m2a = 1.0 - 2*a;

   // save value in double precision
   original_a = a;
   compute_sd();
}

// convert each value to log base sqrt(2)
void krdalign::set_q( blitz::TinyVector<double,4> &newq ) {
   using namespace std;
   for( int i = 0; i < 4; i++ ) {
      q(i) = newq(i);
   }
}

// convert each value to log base sqrt(2)
void krdalign::set_p( blitz::Array<double,2> &newp ) {
   using namespace std;
   assert( newp.rows() == 4 );
   assert( newp.cols() == 4 );

   for( int i = 0; i < 4; i++ ) {
   for( int j = 0; j < 4; j++ ) {
      p(i,j) = newp(i,j);
   }
   }
}

void krdalign::compute_sd() {
   // we want at least this accuracy
   double const ERROR = 0.00001;

   double bottom = 1.0+3.0*original_a - original_a*original_b - original_b;
   double pi_m = 0.5*(1.0-original_b)/bottom;
   double pi_i = original_a/bottom;
   double pi_h = 0.5*original_a*(1.0-original_b)/bottom;
   pr_red = pi_m + pi_i + pi_h;

   // make some assertions here
   assert( pi_m >= 0 && pi_m <= 1.0 );
   assert( pi_i >= 0 && pi_i <= 1.0 );
   assert( pi_h >= 0 && pi_h <= 1.0 );

   // this should be close to one
   assert( std::abs(2*pr_red - 1.0) < ERROR );

   // Now compute the probability that we start in M or I conditioned on the fact that
   // we are already in a red or blue area and hidden states transition immediatly to I.
   pr_m = pi_m/pr_red;
   pr_i = (pi_i+pi_h)/pr_red;
   
   assert( pr_m >= 0 && pr_m <= 1.0 );
   assert( pr_i >= 0 && pr_i <= 1.0 );
   assert( std::abs(pr_m + pr_i - 1.0) < ERROR );
}

// Here we compute the sum off all paths in one of the labeled sub models. Gaps are inserted into sequence B. We compute this
// backwards!
//       0  1  2        
//       B  B  B             
// 0  A  X           
// 1  A  4  2        
// 2  A  7  5  3   
// 3  A  X  8  6     
// 4  A  X  X  9     
//
// Matches are therefore from the left and gaps from the top. We compute
void krdalign::hmmdp( dna_sequence_region_data &a, uint ai, uint aj, dna_sequence_region_data &b, uint bi, uint bj, smatrix &S, int at ) {
   int n = aj - ai + 1, m = bj - bi + 1;
 
    // need these for START positions
   ereal p_m = pr_m, p_i = pr_i;
  
   // a must be longer than b
   assert( n >= m );

#ifdef DEBUG_HMM
    std::cerr << "type\t" << at << std::endl;
#endif

   // compute odds
   odds_a(0) = odds_b(0) = 1.0;
   for( int i = 1; i <= n; i++ ) odds_a(i) = odds_a(i-1) * q((int)a[aj-(uint)(i-1)]);
   for( int j = 1; j <= m; j++ ) odds_b(j) = odds_b(j-1) * q((int)b[bj-(uint)(j-1)]);

   // initialize match at origin and first column, remember we are doing this backwords!
   if( at & END ) { 
      M(1,1) = p((int)a[aj],(int)b[bj]);
      I(1,1) = 0.0;
      S(1,1) = M(1,1) / odds_a(1) / odds_b(1); 
                                    // normally can only accept in hidden state, but we allow 
                                    // accepting from state M here instead

      // if we are finishing an alignment, we can end in an insert state. Doing things backwards,
      // this means we can START in an insert state! If we are ending in an insert state, this means 
      // we are coming from a MATCH state. But things are backwards, so the cost a is added when considering
      // acceptance.

      I(1,0) = q((int)a[aj]);
      M(1,0) = 0.0;
      S(1,0) = I(1,0) / odds_a(1);
      for( int i = 2; i < n; i++ ) {   
         I(i,0) = I(i-1,0) * p_b*q((int)a[aj-i+1]);
         S(i,0) = I(i,0) / odds_a(i);
         M(i,0) = 0.0;
      }

   } else { // unless we are at the end of an alignment, we enter from the hidden state (going backwards)
      M(1,1) = p_a*p((int)a[aj],(int)b[bj]);
      S(1,1) = I(1,1) = 0.0;
      
#ifdef DEBUG_HMM
     std::cerr << "M(1,1)\t" << M(1,1).as_double() << p_a.as_double() << std::endl;
#endif
      // for every other case we exit (start) only from match
      for( int i = 1; i < n; i++ ) {   
         I(i,0) = M(i,0) = S(i,0) = 0.0;
      }
   }

   //std::cerr << "1\t1\t" << M(1,1).as_double() << "\t" << I(1,1).as_double() << std::endl;
   //std::cerr << p_c.as_double() << "\t" << p_e.as_double() << "\t" << p_d.as_double() << std::endl;
   //std::cerr << p_1mamc.as_double() << std::endl;
   //std::cerr << M(1,1).as_double() << "\t" << odds_a(1).as_double() << "\t" <<odds_b(1).as_double() << std::endl;
 
   // set next block to null
   M(1,2) = I(1,2) = S(1,2) = 0.0;
  
   // fill rest of matrix
   for( int j, j_max, i = 2; i <= n; i++ ) {
      j_max = (i <= (int)k) ?  std::min(i,m) : std::min(shape_fx(i-k-1),m);

      for( j = 1; j <= j_max; j++ ) {

        // normally, we end at the other hidden state
         M(i,j) = (M(i-1,j-1)*p_1m2a + I(i-1,j-1)*p_a)*p((int)a[aj-(uint)(i-1)],(int)b[bj-(uint)(j-1)]);
         I(i,j) = (p_b*I(i-1,j) + p_1mb*M(i-1,j)) * q((int)a[aj-(uint)(i-1)]);
         S(i,j) = I(i,j)/odds_a(i)/odds_b(j);

         // If we are starting an alignment, we can end at either of the blue states. We initialize with
         // an adjusted stationary distribution.
         if( (at & START) && (i == n && j == m) ) { 
#ifdef DEBUG_HMM
            std::cerr << "AT START" << std::endl;
#endif
            S(i,j) = (p_m*M(i,j) + p_i*I(i,j))/odds_a(i)/odds_b(j);
         }
#ifdef DEBUG_HMM
         std::cerr << i << "\t" << j;
         std::cerr << "\t" << M(i,j).as_double();
         std::cerr << "\t" << I(i,j).as_double();
         std::cerr << "\t" << S(i,j).as_double();
         std::cerr << std::endl;
#endif
         // set block below to null
         M(i+1,j) = 0.0;
         I(i+1,j) = 0.0;
         S(i+1,j) = 0.0;
      }
      // set next block over to null
      //if( j <= m ) {
        M(i,j) = 0.0;
         I(i,j) = 0.0;
         S(i,j) = 0.0;
      //}

   }
   //std::cerr << "dalign done\n";
   // end should be n,2
}

void krdalign::maindp_init_borders( dparray &dp, dna_sequence_region_data &a, dna_sequence_region_data &b ) {
   int n = (int) a.size(), m = (int)b.size();

   // Top corner is the only valid starting point. Since blue chains to red and red to blue, 
   // we initialize blue to be the probability that we start in red, and red to the probability
   // we start in blue. Since things are symmetric, these are the same.
   init_location(dp,0,0);
   dp(0,0).R.s = dp(0,0).D.s = pr_red;

   // Initialize borders.
   for( int i = 1; i <= n; i++ ) init_location(dp,i,0);
   for( int i = 1; i <= m; i++ ) init_location(dp,0,i);
}

void krdalign::maindp_filldp( dparray &dp, dna_sequence_region_data &adata, dna_sequence_region_data &bdata ) {
   int n = (int) adata.size(), m = (int)bdata.size();
   ereal block_score, best_score, score, match_cost, open_cost;
   for( int i = 1; i <= n; i++ ) {
   for( int j = 1; j <= m; j++ ) {
      init_location(dp,i,j);
      //std::cerr << "ij " << i << "\t" << j << std::endl;
#if DEBUG_DP>0
      if( i == 4 && j == 1 ) {
      std::cerr << "ij " << i << "\t" << j << std::endl;
      std::cerr <<  dp(i,j).D.s.as_double() << std::endl;
      std::cerr <<  dp(i,j).R.s.as_double() << std::endl;
      }
#endif

      // precompute all values for the box
      int a_start = std::max(i-(int)k-(int)shape_z+2,1), b_start=std::max(j-(int)k-(int)shape_z+2,1);
      int at = BETWEEN_BLOCKS;
      if( (a_start == 1) && (b_start == 1) ) at |= START;
      if( (i == n) && (j == m) ) at |= END;
      if( (i-a_start) >= (j-b_start) ) {
         //std::cerr << "dalign " << a_start << " " << i << " " << b_start << " " << j << std::endl;
         if( j-b_start+1 > (int)k ) hmmdp( adata, a_start-1, i-1, bdata, j-k, j-1, Sd, at );
         else hmmdp( adata, a_start-1, i-1, bdata, b_start-1, j-1, Sd, at );

         if( (j-b_start) + 1 > (int)k ) hmmdp( bdata, b_start-1, j-1, adata, i-k, i-1, Sr, at );
         else hmmdp( bdata, b_start-1, j-1, adata, i-(j-b_start)-1, i-1, Sr, at );
      } else {
         //std::cerr << "ralign " << a_start << " " << i << " " << j-(i-a_start) << " " << j << std::endl;
         if( i-a_start + 1 > (int)k ) hmmdp( adata, a_start-1, i-1, bdata, j-k, j-1, Sd, at );
         else hmmdp( adata, a_start-1, i-1, bdata, j-(i-a_start)-1, j-1, Sd, at );

         if( i - a_start + 1 > (int)k )  hmmdp( bdata, b_start-1, j-1, adata, i-k, i-1, Sr, at );
         else hmmdp( bdata, b_start-1, j-1, adata, a_start-1, i-1, Sr, at );
      }

      // look in this box
      for( int a = i; a >= a_start; a-- ) {
      for( int b = j; b >= b_start; b-- ) {
#if DEBUG_DP>0
      if( i == 4 && j == 1 ) std::cerr << "ab " << a << "\t" << b << std::endl;
#endif
         // If we are above the first diagonal we try adding a down block.
         if( (b - j + i) >= a ) {
            if( (i-a+1 > (int)k)  && ((j-b+1) > shape_fx(i-a-k)) ) {
               //std::cerr << "D " <<  a << "\t" << b << std::endl;
               break; // off curve 
            }
            // compute score of down block
            block_score = Sd(i-a+1,j-b+1);
#if DEBUG_DP>1
      if( i == 4 && j == 1 ) {
            std::cerr <<  "D\t" << dp(a-1,b-1).R.s.as_double() << std::endl;
            std::cerr <<  "D\t" << dp(i,j).D.s.as_double() << std::endl;
            std::cerr <<  "D\t" << block_score.as_double() << std::endl;
      }
#endif
            // find optimal chaining, we can only chain from an R state
            score = block_score*dp(a-1,b-1).R.s;
            //std::cerr <<  "D\t" << score.as_double() << std::endl;
            if( score > dp(i,j).D.s ) {
               dp(i,j).D.s = score;
               dp(i,j).D.u = i - a + 1; dp(i,j).D.l = j - b + 1; dp(i,j).D.from = R_block;

            }
               
            // allow chaining of k down blocks gap
            score = block_score*dp(a-1,b-1).D.s;
            if( ((i-a + 1) == (int)k) && ( score > dp(i,j).D.s) ) {
               dp(i,j).D.s = score;
               dp(i,j).D.u = i - a + 1; dp(i,j).D.l = j - b + 1; dp(i,j).D.from = D_block;
            } // if
            
         // Otherwise consider chaining a right block.
         } else {
            // compute score of right block
            if( (i-a+1 > (int)k) || ((j-b+1) > ((int)k+shape_fy(i-a+1))) ) {
               //std::cerr << "R " << a << "\t" << b << std::endl;
               break; // off curve 
            }
            block_score = Sr(j-b+1,i-a+1);
#if DEBUG_DP>1
      if( i == 4 && j == 1 ) {
            std::cerr <<  "R\t" << dp(a-1,b-1).D.s.as_double() << std::endl;
            std::cerr <<  "R\t" << dp(i,j).D.s.as_double() << std::endl;
            std::cerr <<  "R\t" << block_score.as_double() << std::endl;
      }
#endif  
            // find optimal chaining, we can only chain from a D state
            score = block_score*dp(a-1,b-1).D.s;
            if( score > dp(i,j).R.s ) {
               dp(i,j).R.s = score;
               dp(i,j).R.u = i - a + 1; dp(i,j).R.l = j - b + 1; dp(i,j).R.from = D_block;
            }
               
            // allow chaining of k down blocks gap
            score = block_score*dp(a-1,b-1).R.s;
            if( ((j-b + 1) == (int)k) && (score > dp(i,j).R.s) ) {
               dp(i,j).R.s = score;
               dp(i,j).R.u = i - a + 1; dp(i,j).R.l = j - b + 1; dp(i,j).R.from = R_block;
            } // if
         } // if
      } // for b
      } // for a

      // check pure gap end cases
      if( (i == n) && (j == m ) ) {
         ereal gap = 1.0;
         for( int k = 1; k < n; k++ ) {
            ereal score = gap*dp(n-k,m).R.s;
            if( score > dp(n,m).D.s ) {
#if DEBUG_DP>0
               std::cerr << "GAP A\t" << score.as_double() << "\t" << k << std::endl;
#endif
               dp(n,m).D.s = score;
               dp(n,m).D.u = k; dp(n,m).D.l = 0; dp(n,m).D.from = Ga_state;
            }
            gap *= p_b;
         }
            
         gap = 1.0;
         for( int k = 1; k < m; k++ ) {
            ereal score = gap*dp(n,m-k).D.s;
            if( score > dp(n,m).R.s ) {
#if DEBUG_DP>0
               std::cerr << "GAP B\t" << score.as_double() << "\t" << k << std::endl;
#endif
               dp(n,m).R.s = score;
               dp(n,m).R.u = 0; dp(n,m).R.l = k; dp(n,m).R.from = Gb_state;
            }
            gap *= p_b;
         }
      }

      //std::cout << i << '\t' << j << '\t' << (i - dp(i,j).Ru) << '\t' << (j - dp(i,j).Rl) << '\t' << 0 << std::endl;
      //std::cout << i << '\t' << j << '\t' << (i - dp(i,j).Du) << '\t' << (j - dp(i,j).Dl) << '\t' << std::endl;
   } // for j
   } // for i
}

std::pair<similarity_region_deque_ptr,double> krdalign::maindp_traceback( dparray &dp, dna_sequence_region_data &a, dna_sequence_region_data &b ) {
   using namespace std;
   int n = (int) a.size(), m = (int)b.size();

   similarity_region_deque_ptr result = new_similarity_region_deque();
   
   // determine starting block
   int i = n, j = m, next_i = i, next_j = j;
   uint8_t type = ( dp(i,j).R.s < dp(i,j).D.s) ? D_block : R_block, next_type = no_state;
   ereal chain_score = state_score( dp, i, j, type );

   while( (i > 0) || (j > 0) ) {
#ifdef DEBUG_TRACEBACK
      cerr << i << "\t" << j << endl;
#endif
      assert( type != no_state );
      assert( type < num_states );

      switch(type) {
         case Ga_state:
         case R_block:
#ifdef DEBUG_TRACEBACK
            cerr << "R\t" << dp(i,j).R.u << "\t" << dp(i,j).R.l << endl;
#endif
            next_i = i - dp(i,j).R.u;
            next_j = j - dp(i,j).R.l;
            next_type = dp(i,j).R.from;
            break;
         case Gb_state:
         case D_block:
#ifdef DEBUG_TRACEBACK
            cerr << "D\t" << dp(i,j).D.u << "\t" << dp(i,j).D.l << endl;
#endif
            next_i = i - dp(i,j).D.u;
            next_j = j - dp(i,j).D.l;
            next_type = dp(i,j).D.from;
      }
     
      // This is a bit confusing. The Ga and Gb states represent the case where we end in a pure
      // gap. This is possible since we can stop in a gap state. Otherwise, we have at least one
      // match.
      if( next_type == Ga_state ) {
         result->push_front( std::make_pair( closed_interval(next_i,i-1), closed_interval(-1,-1) ) );
      } else if( next_type == Gb_state ) {
         result->push_front( std::make_pair( closed_interval(-1,-1), closed_interval(next_j,j-1) ) );
      } else {
         result->push_front( std::make_pair( closed_interval(next_i,i-1), closed_interval(next_j,j-1) ) );
      }

      i = next_i;
      j = next_j;
      type = next_type;
   }
   return std::make_pair( result, chain_score.as_base() );
}

std::pair<similarity_region_deque_ptr,double> krdalign::maindp( dna_sequence_region_data &a, dna_sequence_region_data &b ) {
   int n = (int) a.size(), m = (int)b.size();
   using namespace std;
   dparray dp(n+1,m+1);
   maindp_init_borders( dp, a, b );
   maindp_filldp( dp, a, b );
   return maindp_traceback( dp, a, b );
}

// default k to 6,give original a and b some values so our asserts do go crazy
krdalign::krdalign() : k(20), shape_y(7), original_a(0.2), original_b(0.8) {
   p.resize(dna_alpha::alphabet.size(),dna_alpha::alphabet.size());
   init_shape();
   init_storage();
   init_pr();
}

void krdalign::set_k( unsigned int new_k ) {
   k = new_k;
   init_shape();
   init_storage();
}

std::pair<similarity_region_deque_ptr,double> krdalign::align( dna_sequence_region &seqa, dna_sequence_region &seqb ) {
   return maindp( seqa.data, seqb.data ); 
}

