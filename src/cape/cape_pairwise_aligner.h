#ifndef CAPE_PAIRWISE_ALIGNER_H__
#define CAPE_PAIRWISE_ALIGNER_H__
#include <boost/tuple/tuple.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <iostream>

//#include "krda_improve.h"
#include "needleman.h"
#include "random_seed_factory.h"
#include "nw_model_parameters.h"
#include "scored_fragment.h"
#include "zero_order_background.h"
#include "sequence_statistics.h"
#include "mm_alignment.h"
#include "mm_trainer.h"
#include "input_data.h"
#include "simple_nw_model_parameters.h"

class cape_pairwise_aligner {
private:
	typedef boost::mt19937 base_generator;
	typedef boost::shared_ptr<base_generator> base_generator_ptr;

	typedef std::vector< std::pair<dna_sequence_region_ptr,dna_sequence_region_ptr> > region_pair_vector;
	typedef boost::shared_ptr<region_pair_vector> region_pair_vector_ptr;

	mm_alignment mma;
	mm_trainer mmt;

	// random number generators
	base_generator_ptr baseg;

	double random_sample;
	unsigned int double_random_sample;

	int extension_algorithm;

	// maximum area we can use the detailed aligner
	double MAX_AREA;

	// filter threshold
	ereal THRESHOLD;

	// threshold at which to stop scoring a fragment
	ereal SHORT_CIRCUIT_THRESHOLD;

	//  max distance from chain, if zero, filter is disabled
	int MAX_DISTANCE_FROM_CHAIN;

	// segmentation extension
	unsigned int EXT;

	// the minimum sized window we use to measure composition
	unsigned int min_composition_window;

	// the number of training passes
	unsigned int training_passes;

	// do we really want to be chatty
	bool verbose;

	// have custom parameters been set
	bool custom_parameters;

	// use custom parameters for anchor scoring or dynamic scoring
	bool custom_scoring;

	// disable the random model?
	bool disable_random_model;

	// number of threads to use
	int nthreads;

   // Background frequency for using a fixed background.
   blitz::TinyVector<double,4> q;

	std::vector<simple_nw_model_parameters> parameter_sets, all_parameter_sets;

	void create_parameter_sets();

	// weight, length, and expected random hits
	typedef boost::tuple<unsigned int, unsigned int, unsigned int> seed_parameters;

	seed_parameters computeOSP( input_data &in, fragment &f );

	scored_fragment_ptr_list_ptr find_anchors( input_data &in, scored_fragment a, bool use_mask = true );
	region_pair_vector_ptr find_segments( input_data &in );

	fragment ensure_minimum_for_composition( input_data &in, fragment &f );

//	triplet_table<logint> compute_conditionals( input_data &in, scored_fragment &f );

	scored_fragment_ptr_vector_ptr chain_filter( scored_fragment_ptr_vector &fragments );
	void score_fragments( input_data &in, scored_fragment_ptr_vector_ptr fragments, bool use_mask );

public:
	cape_pairwise_aligner( int k );	// constructor setting the number of models
	enum ext_alg { FW_LOCAL = 0, SW_LOCAL};

	void set_extension_limit( int limit );
	void set_random_sample( double s );
	void set_required_region_size( unsigned int i );
	void set_active_parameter_sets( unsigned int nsets );
	void set_threshold( double t );
	void set_training_passes( unsigned int t );
	void set_chain_filter_distance( int t );
	void set_num_threads( int n );
	void set_custom_scoring( bool b );
   void set_background(  blitz::TinyVector<double,4> nq );

	pairwise_dna_alignment align( dna_sequence_ptr seqa, dna_sequence_ptr seqb );
	mm_trainer const &train( dna_sequence_ptr seqa, dna_sequence_ptr seqb );
	std::pair<pairwise_dna_alignment,annotation_ptr> align_and_annotate( dna_sequence_ptr seqa, dna_sequence_ptr seqb );
	void write_mma_parameters( std::ostream &out );
	void set_parameters_from_file( std::string filename );
   void initialize_parameters_from_file( std::string filename );
	void mma_disable_random_model();
};

inline pairwise_dna_alignment cape_pairwise_aligner::align( dna_sequence_ptr seqa, dna_sequence_ptr seqb ) {
	std::pair<pairwise_dna_alignment,annotation_ptr> r = align_and_annotate( seqa, seqb );
	return r.first;
}

#endif
