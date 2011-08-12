#include "cape_pairwise_aligner.h"

#include <blitz/tinyvec-et.h>
#include <algorithm>
#include <ext/functional>
#include "parameter_functions.h"
#include "pairwise_dna_alignment.h"
#include "statistics.h"
#include "hash_factory.h"
#include "hash_hit_data.h"
#include "fragment_algorithm.h"
#include "fragment_functions.h"
#include "fragment_selectors.h"
#include "algorithm_extra.h"
#include "alignment_functions.h"
#include "sequence_statistics.h"
#include "sequence_statistics_functions.h"
#include "entropy.h"
#include "io.h"
#include "dynamic_fragment_score_worker.h"
#include "custom_fragment_score_worker.h"

fragment cape_pairwise_aligner::ensure_minimum_for_composition( input_data &in, fragment &f ) {
   fragment new_f = f;
   expand_closed_interval( in.a.seq->data, new_f[0], min_composition_window );
   expand_closed_interval( in.b.seq->data, new_f[1], min_composition_window );
   return new_f;
}

/*
triplet_table<logint> cape_pairwise_aligner::compute_conditionals( input_data &in, scored_fragment &f ) {
   triplet_table<unsigned int> freq;
   if( f[0].a > 0 ) {
      freq = in.a.stats->tri[f[0].b] - in.a.stats->tri[f[0].a-1];
   } else {
      freq = in.a.stats->tri[f[0].b];
   }

   if( f[1].a > 0 ) {
      freq += in.b.stats->tri[f[1].b] - in.b.stats->tri[f[1].a-1];
   } else {
      freq += in.b.stats->tri[f[1].b];
   }

   triplet_table<logint> result;

   for( dna_alpha::alphabet_type::const_iterator i = dna_alpha::alphabet.begin(); i != dna_alpha::alphabet.end(); ++i ) {
   for( dna_alpha::alphabet_type::const_iterator j = dna_alpha::alphabet.begin(); j != dna_alpha::alphabet.end(); ++j ) {
      double asize = 0.0;
      for( dna_alpha::alphabet_type::const_iterator k = dna_alpha::alphabet.begin(); k != dna_alpha::alphabet.end(); ++k ) {
         asize += freq(k->index(),j->index(),i->index());
      }
      for( dna_alpha::alphabet_type::const_iterator k = dna_alpha::alphabet.begin(); k != dna_alpha::alphabet.end(); ++k ) {
         //std::cerr << dna_alpha::decode(*k) << dna_alpha::decode(*j) << dna_alpha::decode(*i) << "\t";
         result(k->index(),j->index(),i->index()) = ((double)freq(k->index(),j->index(),i->index()))/asize;
      }
   }
   }
   return result;
}
*/

void cape_pairwise_aligner::set_active_parameter_sets( unsigned int nsets ) {
	// This list is 1 based. Remember to subtract one for each parameter.
#ifndef VITERBI_EXTENSIONS
	static unsigned int const greedy_opt[] = {106, 88, 150, 62, 59, 77, 103, 145, 44, 121, 132, 31, 84, 93, 72,
		131, 146, 75, 28, 54, 61, 78, 134, 138, 142, 46, 58, 63, 115, 119 };
#else
	// If we are using Viterbi segmentation, then we use a list optimized for SW.
	static unsigned int const greedy_opt[] = {52, 82, 8, 97, 20, 36, 50, 68, 113, 81, 83, 9, 66, 19, 24, 109,
		13, 49, 67, 100, 145, 18, 12, 15, 21, 25, 41, 65, 72, 84};
#endif
	assert( nsets > 0 );
	assert( nsets <= 30 );
	parameter_sets.clear();
	for( unsigned int i = 0; i < nsets; ++i ) {
		parameter_sets.push_back( all_parameter_sets[greedy_opt[i]-1] );
	}
}

void cape_pairwise_aligner::create_parameter_sets() {
	using namespace std;
	vector<double> gap_lengths, gap_opens, distances;
	gap_opens.push_back( 0.03125 ); // -5
	gap_opens.push_back( 0.0625  ); // -4
	gap_opens.push_back( 0.125   ); // -3
	gap_opens.push_back( 0.25    ); // -2

	gap_lengths.push_back( 1.333    ); // -2
	gap_lengths.push_back( 2.0      ); // -1
	gap_lengths.push_back( 3.414214 ); // -0.5
	gap_lengths.push_back( 6.285214 ); // -0.25

	distances.push_back( 0.1 );
	distances.push_back( 0.2 );
	distances.push_back( 0.3 );
	distances.push_back( 0.4 );
	distances.push_back( 0.5 );
	distances.push_back( 0.6 );
	distances.push_back( 0.7 );
	distances.push_back( 0.8 );
	distances.push_back( 0.9 );
	distances.push_back( 1.0 );

	all_parameter_sets.clear();
	for( vector<double>::const_iterator d = distances.begin(); d != distances.end(); ++d )
		for( vector<double>::const_iterator a = gap_opens.begin(); a != gap_opens.end(); ++a )
			for( vector<double>::const_iterator L = gap_lengths.begin(); L != gap_lengths.end(); ++L )
				all_parameter_sets.push_back( simple_nw_model_parameters(*d,*a,*L) );
}

cape_pairwise_aligner::seed_parameters cape_pairwise_aligner::computeOSP( input_data &in, fragment &f  ) {
   using namespace std;
   // expand fragment to minimum size
   fragment fcomp = ensure_minimum_for_composition( in, f );

   // estimate region composition
   blitz::TinyVector<double,4> q = composition_frequency(*(in.a.stats),fcomp[0]) + composition_frequency(*(in.b.stats),fcomp[1]);
   q /= (double) (fcomp[0].size() + fcomp[1].size() );

   // compute match probability for a column
   double pr_match = dot(q,q);
   double area = (double)f[0].size() * (double)f[1].size();
   unsigned int w = (unsigned int) (log(random_sample/area)/log(pr_match));
   unsigned int l = (unsigned int) ((double)w*1.4);
   if( l < 5 ) l = 5;      // minimum length for random seeds is 5
   if( w < 3 ) w = 3;      // minimum weight 3
   unsigned int expected_random_hits = (unsigned int)(pow( pr_match, (double)w )*area);
   return seed_parameters(w,l,expected_random_hits);
}

cape_pairwise_aligner::cape_pairwise_aligner( int k ) :	mma(k), mmt(k),
random_sample(2000),
double_random_sample(4000),
extension_algorithm(cape_pairwise_aligner::FW_LOCAL),
MAX_AREA(500*500),
MAX_DISTANCE_FROM_CHAIN(200),
EXT(200),
min_composition_window(2000),
training_passes(5),
verbose(false),
custom_scoring(false),
disable_random_model(false),
nthreads(1)
{
   q = 0.25, 0.25, 0.25, 025;
	THRESHOLD.set_base(70.0);
	SHORT_CIRCUIT_THRESHOLD.set_base(200.0);

	create_parameter_sets();
	set_active_parameter_sets( 10 );
	entropy random;
   unsigned int rseed = (unsigned int) random.get_int();
   if( rseed == 0 ) rseed = 1; // bug in old boost doesn't allow 0 as a seed
   base_generator::result_type brseed(rseed);
   baseg = base_generator_ptr( new base_generator( brseed ) );
}


void cape_pairwise_aligner::set_required_region_size( unsigned int i ) {
   MAX_AREA = (double)i*(double)i;
}

void cape_pairwise_aligner::set_random_sample( double rs ) {
   random_sample = rs;
	double_random_sample = (unsigned int)rs*2;
}


void cape_pairwise_aligner::set_custom_scoring( bool b ) { 
   custom_scoring = b;
}

void cape_pairwise_aligner::set_threshold( double t ) { 
   THRESHOLD.set_base(t);
}

void cape_pairwise_aligner::set_training_passes( unsigned int t ) {
	training_passes = t;
}

void cape_pairwise_aligner::set_chain_filter_distance( int t ) {
	MAX_DISTANCE_FROM_CHAIN = t;
}

void cape_pairwise_aligner::set_num_threads( int n ) {
	nthreads = n;
}

void cape_pairwise_aligner::set_extension_limit( int limit ) {
   //fal.set_extension_limit(limit);
   //swl.set_extension_limit(limit);
}

mm_trainer const &cape_pairwise_aligner::train( dna_sequence_ptr seqa, dna_sequence_ptr seqb ) {
	using namespace std;

	// preprocess and package input data
   input_data indata( seqa, seqb );

	region_pair_vector_ptr segments = find_segments( indata );

	if( disable_random_model ) mmt.disable_random_model();

	cerr << "INITIAL PARAMETERS" << endl << endl;
	mmt.pretty_print_parameters( cerr );

	// now we train on all the regions
	cerr  << endl << "Training on " << segments->size() << " regions." << endl;

	ereal prev_sc = 0.0;
	for( unsigned int j = 0; j < training_passes; ++j ) {
		cerr << endl << "Pass " << (j+1) << " of " << training_passes << "." << endl;
		mmt.start_training();
		mmt.constrained_bw_train( segments->begin(), segments->end(), nthreads );
		ereal pr = mmt.end_training();

		cerr << "Forward pr:" << pr.as_base() << ", " << (double)pr << endl;
		mmt.write_parameters( cerr );

		// tie parameters
		for( int k = 0; k < mmt.num_models(); ++k ) {
			nw_model_parameters_ptr nwp = mmt.get_parameters(k);
			mmt.set_parameters(k,*nwp,(int)mmt.get_region_length(k));
		}
		mmt.set_random_param((int)mmt.get_random_length_a(), (int)mmt.get_random_length_b(), mmt.get_random_distribution());

		// check if we are finished
		if( fabs( (pr.as_base()-prev_sc.as_base()) ) < 1 ) break;
		prev_sc = pr;
	}

	return mmt;
}

cape_pairwise_aligner::region_pair_vector_ptr cape_pairwise_aligner::find_segments( input_data &indata ) {
	using namespace std;
	// start with one large region (implicit anchors at each end)
   scored_fragment_ptr_list anchors;
   scored_fragment_ptr front(new scored_fragment()), back(new scored_fragment());
	// we use 0,0 so that the midpoint is a valid sequence position
   front->add( closed_interval(0,0) );
   front->add( closed_interval(0,0) );
	// right side is always taken to be midpoint-1, so this indicates the end of the sequence
   back->add( closed_interval(indata.a.seq->data.size(),indata.a.seq->data.size()) );
   back->add( closed_interval(indata.b.seq->data.size(),indata.b.seq->data.size()) );
   anchors.push_back( front );
   anchors.push_back( back );


   double random_sample_original = random_sample;
   scored_fragment_ptr_list_iter prev = anchors.begin(), i = ++anchors.begin();
   while ( i != anchors.end() ) {
      // if area too small, get new anchors for it
      if( interior_area( **prev, **i ) > MAX_AREA ) {
         // if either sequence has >=40% repeats, then don't mask them out
         scored_fragment ifrag = interior( **prev, **i );
         //bool use_mask =
			//(fraction_repeats(*(indata.a.stats), ifrag[0]) >= 0.40  || fraction_repeats(*(indata.b.stats), ifrag[1]) >= 0.40)
			//? false : true;
         scored_fragment_ptr_list_ptr new_anchors = find_anchors( indata, interior( **prev, **i ), true );
         anchors.splice( i, *new_anchors );

         // adjust i
         i = prev;
         ++i;

			// small enough, skip to next range
      } else {
         random_sample = random_sample_original;
         ++i;
         ++prev;
      }
   }
	// now filter out as many guessed fragments as possible
	if( anchors.size() > 2 ) {
		prev = anchors.begin();
		i = ++anchors.begin();
		scored_fragment_ptr_list_iter j = i;
		++j;

		while( j != anchors.end() ) {
			if( (**i).score() < THRESHOLD && interior_area(**prev,**j) <= MAX_AREA ) {
				i = anchors.erase(i);
				++j;
			} else {
				++prev;
				++i;
				++j;
			}
		}
	}

   cerr << endl << "Finished anchoring with " << anchors.size() << " anchors." << endl << endl;

	// first we create a set of sequence regions
	region_pair_vector_ptr segments( new region_pair_vector() );

   for( prev = anchors.begin(), i = ++anchors.begin(); i != anchors.end(); ++i, ++prev ) {
      scored_fragment region = midpoint_interior( **prev, **i );
		if( verbose ) cerr << region << endl;
      assert( region[0].a < region[0].b );
      assert( region[1].a < region[1].b );

		dna_sequence_region_ptr region_a( new dna_sequence_region(indata.a.seq,region[0]) );
		dna_sequence_region_ptr region_b( new dna_sequence_region(indata.b.seq,region[1]) );
		segments->push_back( make_pair(region_a,region_b) );
	}
	return segments;

}

std::pair<pairwise_dna_alignment,annotation_ptr> cape_pairwise_aligner::align_and_annotate( dna_sequence_ptr seqa, dna_sequence_ptr seqb ) {
	using namespace std;

	// preprocess and package input data
   input_data indata( seqa, seqb );

	region_pair_vector_ptr segments = find_segments( indata );

	if( disable_random_model ) mma.disable_random_model();

	// do alignment
	cerr << "Aligning " << segments->size() << " regions." << endl;

	annotation_ptr ann;
	dna_alignment_sequence_ptr alignmenta( new dna_alignment_sequence() );
   dna_alignment_sequence_ptr alignmentb(new dna_alignment_sequence() );
   pairwise_dna_alignment final = pairwise_dna_alignment( alignmenta, alignmentb, 0.0 );

	mma.set_interface_size(4);
	tie(final,ann) = mma.constrained_align_and_annotate( segments->begin(), segments->end() );
	cerr << "Done." << endl;
   return make_pair( final, ann );
}


// here we filter by distance to optimal chain
scored_fragment_ptr_vector_ptr cape_pairwise_aligner::chain_filter( scored_fragment_ptr_vector &fragments ) {
	using namespace std;

	// perhaps the problem with sorting is that the reference counting gets borked and thus a fragment
	// is accidentally deleted. If so, keeping a copy of each pointer here should stop it.
	scored_fragment_ptr_vector copy_test = fragments;

	// hold the final result
	scored_fragment_ptr_vector_ptr result( new scored_fragment_ptr_vector() );
	assert( fragments.size() > 0 );

	ereal basic_score;
	basic_score.set_base(200.0);
	for( scored_fragment_ptr_vector_iter i = fragments.begin(); i != fragments.end(); ++i ) (**i).set_score(basic_score );

   // first find optimal chain
	fragment_algorithm::optimal_chain( fragments, back_inserter(*result), fragment_algorithm::gap_score_cf(1.0 - 1.0/mma.get_random_length_a()) );
	//fragment_algorithm::optimal_chain_fast( fragments, back_inserter(*result) );
//	fragment_algorithm::optimal_chain( fragments, back_inserter(chain) );

	//for( scored_fragment_ptr_vector_iter i = fragments.begin(); i != fragments.end(); ++i ) assert( i->get() != 0 );
	assert( result->size() > 0 );
	return result;
}

void cape_pairwise_aligner::score_fragments( input_data &in, scored_fragment_ptr_vector_ptr fragments, bool use_mask ) {
	using namespace std;

	// score each fragment
	concurrent_queue<scored_fragment_ptr> fqueue;
	copy( fragments->begin(), fragments->end(), back_inserter(fqueue) );
	fqueue.shutdown_when_empty();
	boost::thread_group workers;
	
	vector<nw_model_parameters> nwsets;
	for( int i = 0; i < mma.num_models(); ++i ) nwsets.push_back( mma.get_parameters(i) );
	custom_fragment_score_worker worker( nwsets, in, fqueue );
	for( int i = 0; i < nthreads; ++i ) workers.create_thread(worker);

   workers.join_all();
}


scored_fragment_ptr_list_ptr cape_pairwise_aligner::find_anchors( input_data &in, scored_fragment f, bool use_mask ) {
    using namespace std;
    using namespace __gnu_cxx;
    using namespace boost;
    
    seed_parameters p = computeOSP(in,f);
    random_seed_factory sf( get<1>(p) );
    seed myseed = sf.create( get<0>(p) );

    hash_factory<hash_hit_data,dna_sequence> hsh( myseed );
    hsh.set_use_mask( use_mask );
    hash_hit_data_ptr sequence_hash = hsh.create( in.a.seq, in.b.seq, f );
    scored_fragment_ptr_vector_ptr all_fragments( new scored_fragment_ptr_vector() );
    fragment_algorithm::create_fragments( *sequence_hash, back_inserter(*all_fragments) );
    scored_fragment_ptr_vector_ptr fragments = all_fragments;

    // score each fragment
    cerr << f << " w: " << get<0>(p) << " l: " << get<1>(p) << " expected/found: " << get<2>(p) << "/" << fragments->size()
		  << " p: " << myseed.pattern_string();
    if( use_mask ) cerr << " repeats masked,";
    cerr << " scoring..." << endl;

    if( all_fragments->size() == 0 ) throw runtime_error("cape_pairwise_anchor: found no fragments");

    // use chain filter to reduce number of test points
    if( MAX_DISTANCE_FROM_CHAIN > 0 ) {
		fragments = chain_filter( *all_fragments );
		cerr << "filtered fragments: " << (all_fragments->size()-fragments->size()) << endl;
	}

    // clear score on each fragment
	for( scored_fragment_ptr_vector_iter i = all_fragments->begin(); i != all_fragments->end(); ++i ) (**i).set_score(0.0);

	score_fragments( in, fragments, use_mask );

	// split fragments into good ones and poor scoring ones
	scored_fragment_ptr_vector good_fragments;
	scored_fragment_ptr_vector poor_fragments;
	split_copy( *fragments, score_below(THRESHOLD), back_inserter(poor_fragments), back_inserter(good_fragments) );

   // filter out any fragments with score below the sample deviation
   cerr << "before filter: " << fragments->size();
   cerr << "\tafter filter: " << good_fragments.size();

  /* for( scored_fragment_ptr_vector_iter i = good_fragments.begin(); i != good_fragments.end(); ++i ) {
      cerr << **i << "\t" <<(**i).score() << endl;
   }
   for( fragment_ptr_vector_iter i = fragments->begin(); i != fragments->end(); ++i ) {
      cerr << (**i)[0].a << "\t" <<(**i)[1].a << endl;
   }
   exit(1);*/

   // chain fragments
   scored_fragment_ptr_list_ptr chain( new scored_fragment_ptr_list() );
	fragment_algorithm::optimal_chain( good_fragments, front_inserter(*chain), fragment_algorithm::null_cf() );

	cerr << "\tchain size: " << chain->size() << endl;

	// Now we guess if we must. For each region in the chain, if we will end up
	// using the same parameter set as we just did, we may as well just guess now
	// and avoid duplicate work. We keep guessing until either the tiling is dense
	// enough or we will use a new seed on the next pass.

	// First we add a start and end fragment to our chain.
	pair<scored_fragment_ptr,scored_fragment_ptr> ends = build_ends(f);
	chain->push_front(ends.first);
	chain->push_back(ends.second);


	// iterate over each area
	for( scored_fragment_ptr_list_iter i = ++(chain->begin()), prev = chain->begin(); i != chain->end(); ) {
		// skip if area is small enough
		if( interior_area( **prev, **i ) <= MAX_AREA ) { ++prev; ++i; continue; }

		// otherwise get the parameters we will use for it
		fragment ifrag = interior( **prev, **i );
		seed_parameters next_p = computeOSP(in,ifrag);

		// if we will use a new weight, leave for next round
		if( get<0>(next_p) != get<0>(p) ) { ++prev; ++i; continue; }

		// Otherwise we need to choose a random fragment in the region from all_fragments.
		// First we filter out fragments that can't be in the set.
		scored_fragment_ptr_vector_ptr filtered( new scored_fragment_ptr_vector() );
		for( scored_fragment_ptr_vector_iter j = all_fragments->begin(); j != all_fragments->end(); ++j )
			if( (**prev >> **j) && (**j >> **i) ) filtered->push_back( *j );

		// We didn't find any good fragment to insert, and will repeat this search next time.
		// So we are really out of luck now. This should never happen of course, but just in case.
		if( filtered->size() == 0 ) {
			throw runtime_error("cape_pairwise_anchor: found no fragments causing endless loop");
		}

		if( verbose ) cerr << "Rescoring a possible " << filtered->size() << " fragments for filling." << endl;
		// Score any unscored fragments.
		score_fragments(in, filtered, use_mask);

		/*
		// Try to chain again.
		scored_fragment_ptr_list_ptr subchain( new scored_fragment_ptr_list() );
		fragment_algorithm::optimal_chain( *filtered, front_inserter(*subchain), fragment_algorithm::null_cf() );

		if( subchain->size() > 0 ) {
			cerr << "chained: " << subchain->size() << endl;

		    chain->splice( i, *subchain );
			i = prev;
			++i;
		} else {
			// SHOULD NEVER GET HERE
			assert(0 && "error in anchoring" );
*/
			// We sort the fragments by score.
			sort( filtered->begin(), filtered->end(), by_increasing_fragment_score() );
	//		cerr << filtered->back()->score() << "\t" << filtered->front()->score() << endl;

			// choose the best scoring fragment
			if( verbose ) cerr << "greedily chose " << *(filtered->back()) << endl;
			i = chain->insert( i, filtered->back() );
//		}
	}

	// remove the front and back fake anchors
	chain->pop_front(); chain->pop_back();

	// output and return our final chain
   if( verbose ) for( scored_fragment_ptr_list_citer i = chain->begin(); i != chain->end(); ++i ) {
      cerr << **i << "\t" << (**i).parameter_set() << endl;
   }

   return chain;
}


void cape_pairwise_aligner::set_parameters_from_file( std::string filename ) {
   using namespace std;
   ifstream in( filename.c_str() );
   if( in.fail() ) throw runtime_error("set_parameters_from_file: could not open parameter file");

	vector<double> mdist;
	mdist.push_back( read<double>(in) );
	int rand_a = (int)read<double>(in);
	int rand_b = (int)read<double>(in);
	blitz::TinyVector<double,4> q;
	q = read<double>(in), read<double>(in), read<double>(in), read<double>(in);
	mma.set_random_param(rand_a, rand_b, q);
	mmt.set_random_param(rand_a, rand_b, q);

	int nmodels = read<int>(in);
	for( int j = 0; j < nmodels; ++j ) {
		mdist.push_back( read<double>(in) );
		int length = (int)read<double>(in);
		nw_model_parameters param;
		in >> param;
		mma.set_parameters( j, param, length );
		mmt.set_parameters( j, param, length );
	}
	mma.set_region_dist( mdist );
	mmt.set_region_dist( mdist );
}

void cape_pairwise_aligner::initialize_parameters_from_file( std::string filename ) {
   using namespace std;
   ifstream infile(filename.c_str());
	if( infile.fail() ) throw runtime_error("initialize_parameters_from_file: could not open initialization file");
	unsigned int nmodels = read<unsigned int>(infile);

	for( unsigned int k = 0; k < nmodels; ++k ) {
		nw_model_parameters newp;
		double d = read<double>(infile);
		newp.pr_open = read<double>(infile);
		newp.mean_gap_length = read<double>(infile);
		unsigned int explicit_background = read<unsigned int>(infile);

		if( explicit_background ) {
			newp.q(0) = read<double>(infile);
			newp.q(1) = read<double>(infile);
			newp.q(2) = read<double>(infile);
			newp.q(3) = read<double>(infile);
		} else newp.q = q;
		newp.set_p_hky( d, 1.0 );
		mma.set_parameters( k, newp, 200 );
		mmt.set_parameters( k, newp, 200 );
	}

	infile.close();
}

void cape_pairwise_aligner::mma_disable_random_model() {
	disable_random_model = true;
}

void cape_pairwise_aligner::set_background(  blitz::TinyVector<double,4> nq ) {
   for( int i = 0; i < 4; ++i ) q(i) = nq(i);
}

