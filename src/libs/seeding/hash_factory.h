#ifndef HASH_FACTORY_H__
#define HASH_FACTORY_H__

#include <iostream>
#include "seed.h"
#include "sequence_region.h"

template < class HashResult, class Sequence >
class hash_factory {
   protected:
      // Define some type information.
      typedef typename Sequence::alphabet alpha;
      typedef boost::shared_ptr<Sequence> sequence_ptr;
      typedef typename Sequence::data_type sequence_data;
      typedef typename sequence_data::const_iterator sequence_data_citer;
      typedef std::vector<sequence_ptr> sequence_ptr_vector;
      typedef sequence_region<sequence_ptr> my_sequence_region;
      typedef typename my_sequence_region::data_type my_sequence_region_data;
      typedef typename my_sequence_region_data::const_iterator my_sequence_region_data_citer;
      typedef boost::shared_ptr<my_sequence_region> my_sequence_region_ptr;
      typedef std::vector<my_sequence_region_ptr> my_sequence_region_ptr_vector;

      typedef unsigned int uint;

      // Seed we hash with.
      seed s;

      // We keep this around to reduce allocations.
      sequence_data match;

      // use masking information
      bool use_mask;

		// generate a match of a spaced seed to a sequence
     	sequence_data const &seed_match( my_sequence_region_data_citer data ) {
         using namespace std;

         // start with empty sequence
         match.clear();

         for( uint j = 0; j < (uint)s.length(); j++, data++ ) {
            // ignore regions containing repeats
            if( data->masked() && use_mask ) {
               match.clear();
               return match;
            }
            if( s.pattern_string()[j] == '1' )  match.push_back( *data );
         } // for
         return match;
      } // function

		// hash a sequence into the hash data
		void hash_sequence( HashResult &data, uint seq_ind, my_sequence_region_ptr seq_ptr ) {
			using namespace std;
		   my_sequence_region &seq = *seq_ptr;

		   // break if region is too small to use seed
		   if( seq.data.size() < s.length() )  return;

		   // get start sequence data (sequence starts at 0)
		   my_sequence_region_data_citer i = seq.data.begin();
		   int pos = seq.data.range.a;

		   // find last possible match position + 1 (prevent seed from going past end)
		   my_sequence_region_data_citer end = seq.data.begin() + seq.data.size() - s.length() + 1;

		   for( ; i != end; i++, pos++ ) {
            sequence_data const &my_match = seed_match( i );
            if( my_match.size() > 0 ) data.add( seq_ind, pos, my_match );
		   } // for
		}

   public:

		//--------------------------------------------------------------------------
		// Constructors.
		//--------------------------------------------------------------------------

		hash_factory() : use_mask(false) {}

		hash_factory( seed &newSeed ) : s( newSeed ), use_mask(false) { match.reserve(s.length()); }

      void set_use_mask( bool um ) { use_mask = um; }

		boost::shared_ptr<HashResult> create( sequence_ptr seqa, sequence_ptr seqb, fragment &f ) {
         // bundle sequences into sequence regions
         my_sequence_region_ptr_vector seqList;
         my_sequence_region_ptr newa( new my_sequence_region(seqa,f[0]));
         my_sequence_region_ptr newb( new my_sequence_region(seqb,f[1]));
         seqList.push_back( newa );
         seqList.push_back( newb );
         return this->create(seqList);
      }

      // takes a pair of sequences
		boost::shared_ptr<HashResult> create( sequence_ptr seqa, sequence_ptr seqb ) {
         // bundle sequences into sequence regions
         my_sequence_region_ptr_vector seqList;
         my_sequence_region_ptr newa( new my_sequence_region( seqa ) ), newb( new my_sequence_region( seqb ) );
         seqList.push_back( newa );
         seqList.push_back( newb );
         return this->create(seqList);
      }

		boost::shared_ptr<HashResult> create( my_sequence_region_ptr_vector &seqList ) {
		   // create new hash and set match length
		   boost::shared_ptr<HashResult> data(new HashResult( s.length(), seqList.size() ));

		   // hash each sequence
		   for( uint i = 0; i < seqList.size(); i++ ) {
		      this->hash_sequence( *data, i, seqList[i] );
		   }
		   return data;
		}

      inline void set_seed( seed &newSeed ) {
         s = newSeed;
         match.reserve(s.length());
      }
};

#endif
