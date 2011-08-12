/** Local Extender Search Tool (LEAST)
 *
 *  This is a controller for a generic local alignment. The specific behavior
 *  depends on the EXTENDER, ALIGNER, and TRAINER chosen. Ideally, these three would
 *  be a single class, but due to changes in design during development they are
 *  separate.
 */
#ifndef LEAST_H_
#define LEAST_H_
#include "dna_sequence.h"
#include "extension.h"
#include "scored_fragment.h"
#include "input_data.h"
#include "transformed_fragment.h"
#include "random_seed_factory.h"
#include "fast_hash_factory.h"
#include "fast_hash_hit_data.h"
#include "fragment_algorithm.h"
#include "fragment_functions.h"
#include "input_data.h"
#include "point.h"
#include "point_hash.h"
#include "transformed_fragment.h"
#include "reverse_dna_sequence_region.h"
#include "local_chain.h"
#include "chaining_algorithm.h"
#include "point_hit_builder.h"
#include "sequence_statistics.h"
#include "sequence_statistics_functions.h"
#include "utility.h"

#include <blitz/tinyvec.h>
#include <boost/lexical_cast.hpp>
#include <btl/fasta_writer.h>
#include <hash_map>
#include <ctime>

template<class EXTENDER, class ALIGNER, class TRAINER>
class least {
private:
	EXTENDER &extender;
	ALIGNER &aligner;
	TRAINER &trainer;

	/** Target sequence. */
	dna_sequence_ptr target;
	sequence_statistics_ptr target_stats;

	/** Frequencies for target sequence. */
	blitz::TinyVector<unsigned int,4> countT;

	/** X-drop threshold. */
	double drop_t;

	/** Ungapped x-drop threshold. */
	double udrop_t;

	/** Minimum extension score. */
	double ext_cutoff;

	/** Minimum ungapped extension score. */
	double uext_cutoff;

	/** Currently in training mode. */
	bool   training_mode;

	/** Maximum number of threads to use. */
	unsigned int num_threads;

	/** We output results here. */
	std::ostream *output_stream;

	/** Use viterbi cut. */
	bool use_viterbi_cut;

	/** Viterbi cut threshold. */
	double viterbi_cut_threshold;

	/** Use viterbi training. */
	bool use_viterbi_training;

	/** Find hit points between target and query. */
	boost::shared_ptr< std::vector<point> > find_hit_points( dna_sequence_ptr query );

	/** Align the extension. */
	bool align( dna_sequence_ptr, extension &ext );

	/** Train on the extension. */
	void train( dna_sequence_ptr, extension &ext );

	typedef std::vector<std::pair<dna_sequence_region_ptr,dna_sequence_region_ptr> > region_pair_vector;
	typedef boost::shared_ptr<region_pair_vector> region_pair_vector_ptr;

public:
	least(EXTENDER &ex, ALIGNER &a, TRAINER &mt);

	/** Set the alignment parameters from a file */
	void set_parameters_from_file( std::string filename );

	/** Set the alignment parameters from a file */
	void initialize_parameters_from_file( std::string filename );

	/** Set the drop threshold for extensions. */
	void set_drop_threshold(double);

	/** Set the cut off threshold for extensions. */
	void set_cut_threshold(double);

	/** Set the drop threshold for ungapped extensions. */
	void set_ungapped_drop_threshold(double);

	/** Set the cut off threshold for ungapped extensions. */
	void set_ungapped_cut_threshold(double);

	/** Set stop at repeat behaviour. */
	void set_stop_at_repeat( bool f );

	/** Maximum number of threads to use. */
	void set_num_threads( unsigned int t );

	/** Output stream for results. */
	void set_output_stream( std::ostream *o );

	/** Set viterbi cut threshold. */
	void set_viterbi_cut( double c );

	/** Toggle viterbi training. */
	void set_viterbi_training( bool t );

	/** Perform a local alignment of query against target. */
	void search( dna_sequence_ptr query );

	/** Set target sequence. */
	void set_target( dna_sequence_ptr t );

	/** Start a training pass. */
	void start_training();

	/** Stop training and adjust parameters. */
	ereal stop_training();

	/* Outdated methods.
	void set_hcost(double);
	void set_vcost(double);
	void set_min_extend_score(double);
	void set_min_chain_score(double);
	void set_anchor_search_step(int);
	void set_anchor_search_max_box(int);
	*/
};

template<class EXTENDER, class ALIGNER, class TRAINER>
least<EXTENDER,ALIGNER,TRAINER>::least(EXTENDER &ex, ALIGNER &al, TRAINER &tr ) :
	extender(ex),
	aligner(al),
	trainer(tr),
	drop_t(-25.0),
	udrop_t(-10.0),
	ext_cutoff(10.0),
	uext_cutoff(10.0),
	training_mode(false),
	num_threads(1),
	use_viterbi_cut(false),
	viterbi_cut_threshold(0.0),
	use_viterbi_training(false) {}


template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_drop_threshold( double t ) { drop_t = t; }

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_ungapped_drop_threshold( double t ) { udrop_t = t; }

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_cut_threshold( double t ) { ext_cutoff = t; }

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_ungapped_cut_threshold( double t ) { uext_cutoff = t; }

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_stop_at_repeat( bool t ) { extender.set_stop_at_repeat(t); }

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_num_threads( unsigned int t ) { num_threads = t; }

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_output_stream( std::ostream *o ) { output_stream = o; }

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_parameters_from_file( std::string filename ) {
	using namespace std;

	ifstream infile(filename.c_str());
	if( infile.fail() ) throw runtime_error("set_parameters_from_file: could not open parameter file");
	aligner.read_parameters(infile);
	infile.seekg(0,ios::beg);
	extender.read_parameters(infile);
	infile.seekg(0,ios::beg);
	trainer.read_parameters(infile);
	infile.close();
}

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::initialize_parameters_from_file( std::string filename ) {
	using namespace std;

	ifstream infile(filename.c_str());
	if( infile.fail() ) throw runtime_error("initialize_parameters_from_file: could not open initialization file");
	unsigned int nmodels = read<unsigned int>(infile);

	// Get default background.
	unsigned int total = sum(countT);
	blitz::TinyVector<double,4> q;
	for(int i = 0; i < 4; ++i ) q(i) = (double)countT(i)/(double)(total);

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
		aligner.set_parameters( k, newp, 200 );
		extender.set_parameters( k, newp, 200 );
		trainer.set_parameters( k, newp, 200 );
	}

	infile.close();
}


template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_target( dna_sequence_ptr t ) {
	target = t;
	target_stats = measure_statistics( *target );
	countT = composition_frequency(*target_stats,closed_interval(0,target->data.size()-1));
	unsigned int total = sum(countT);
	blitz::TinyVector<double,4> q;
	for(int i = 0; i < 4; ++i ) q(i) = (double)countT(i)/(double)(total);

	extender.set_default_parameters(q);
	extender.disable_random_model();

	aligner.set_default_parameters(q);
	aligner.disable_random_model();

	trainer.set_default_parameters(q);
	trainer.disable_random_model();
}


template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_viterbi_cut( double t ) {
	viterbi_cut_threshold = t;
}

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_viterbi_training( bool t ) {
	use_viterbi_training = t;
}


template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::start_training() {
	training_mode = true;
	if( use_viterbi_training ) aligner.start_training();
	else trainer.start_training();
}

template<class EXTENDER, class ALIGNER, class TRAINER>
ereal least<EXTENDER,ALIGNER,TRAINER>::stop_training() {
	using namespace std;
	training_mode = false;
	ereal sc = 0.0;
	if( use_viterbi_training ) {
		sc.set_base( aligner.end_training() );
		// Make sure the extender/bw-trainer has the learned parameters set.
		vector<double> dist;
		dist.push_back(0.0);
		for( int k = 0; k < aligner.num_models(); ++k ) {
			extender.set_parameters( k,	aligner.get_parameters(k), aligner.get_region_length(k) );
			dist.push_back(aligner.get_region_pr(k));
		}
		extender.set_region_dist( dist );
	} else {
		sc = trainer.end_training();
	}
	return sc;
}

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::search( dna_sequence_ptr query ) {
	using namespace std;
	using namespace boost;
	using namespace btl;
	using namespace blitz;

	list<extension_ptr> good_extensions;
	extender.reset();

	// Compute combined background.
	sequence_statistics_ptr query_stats = measure_statistics( *query );
	TinyVector<unsigned int,4> countQ = composition_frequency(*query_stats,closed_interval(0,query->data.size()-1));
	TinyVector<unsigned int,4> totalCount = countQ;
	totalCount += countT;
	unsigned int total = sum(totalCount);
	TinyVector<double,4> q;
	for(int i = 0; i < 4; ++i ) q(i) = (double)totalCount(i)/(double)(total);
	zero_order_background oddsT( q, target->data ), oddsQ( q, query->data );

#ifdef TIMING_OUTPUT
	clock_t start_time, end_time;
	start_time = clock();
#endif

#ifndef ANCHOR_AT_ORIGIN
	shared_ptr< vector<point> > hit_points = find_hit_points(query);
#endif

#ifdef TIMING_OUTPUT
	end_time = clock();

	cerr << "seed_time " << (double)(end_time - start_time)/(double)CLOCKS_PER_SEC << endl;
#endif

#ifdef HIT_OUTPUT
#ifndef ANCHOR_AT_ORIGIN
	cerr << "total_hits " << hit_points->size() << endl;
#endif
#endif

	deque<point> seed_points;
	deque<ereal> point_score;

#ifdef TIMING_OUTPUT
	start_time = clock();
#endif

#ifndef ANCHOR_AT_ORIGIN
	// Since the default region lengths won't work well for the ungapped case, we save them
	// and just set everything to 50, since we expect these alignments to be short for difficult
	// cases.
	vector<nw_model_parameters_ptr> mm_params;
	vector<int> mm_lengths;
	for( int k = 0; k < extender.num_models(); ++k ) {
		mm_params.push_back( extender.get_parameters(k) );
		mm_lengths.push_back( (int)extender.get_region_length(k));
		extender.set_parameters(k, *(mm_params.back()), 50 );
	}

	for( vector<point>::iterator i = hit_points->begin(); i != hit_points->end(); ++i ) {
		ereal sc = extender.ungapped_extend_from_point( target, oddsT, query, oddsQ, *i, udrop_t );
		if( sc.as_base() >= uext_cutoff ) {
			seed_points.push_back(*i);
			point_score.push_back(sc);
		}
	}

	for( int k = 0; k < extender.num_models(); ++k )  extender.set_parameters(k, *(mm_params[k]), mm_lengths[k] );
#endif

#ifdef TIMING_OUTPUT
	end_time = clock();
	cerr << "ungapped_filter_time " << (double)(end_time - start_time)/(double)CLOCKS_PER_SEC << endl;
#endif

#ifdef HIT_OUTPUT
	cerr << "filtered_hits " << seed_points.size() << endl;
#endif

#ifdef ANCHOR_AT_ORIGIN
	seed_points.push_back(point(0,0));
	point_score.push_back(0.0);
#endif

/*
	seed_points.clear();
	point_score.clear();
	seed_points.push_back(point(211,99790));
	point_score.push_back(0.0);
*/

#ifdef TIMING_OUTPUT
	start_time = clock();
#endif

	while( ! seed_points.empty() ) {
		point p = seed_points.front();
		seed_points.pop_front();

		ereal cscore = point_score.front();
		point_score.pop_front();

#ifdef DEBUG_OUTPUT
		cerr << "exploring " << p << " " << cscore.as_base() << endl;
#endif

		// do forward extension
		extension_ptr myext(new extension);

		//		mme.extend_from_point(seqA, oddsA, seqB, oddsB, p, back_inserter(extension), drop_t );
		myext->score = extender.extend_from_point(target, oddsT, query, oddsQ, p, back_inserter(myext->segments), drop_t, ext_cutoff );

#ifdef DEBUG_OUTPUT
		for (fragment_ptr_list::iterator i = myext->segments.begin(); i != myext->segments.end(); ++i) {
			cerr << (**i)[0] << "\t" << (**i)[1] << endl;
		}
		cerr << "---" << endl;
#endif

		if( myext->segments.size() > 0 && myext->score.as_base() >= ext_cutoff )
			if( training_mode ) train( query, *myext );
			else align( query, *myext );

	}
#ifdef TIMING_OUTPUT
	end_time = clock();
	cerr << "extend_time " << (double)(end_time - start_time)/(double)CLOCKS_PER_SEC << endl;
#endif


}

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::train( dna_sequence_ptr query, extension &ext ) {
	using namespace boost;
	using namespace std;

	// Skip if we don't pass viterbi cut threshold.
	if( use_viterbi_cut && !align( query, ext) ) return;

	// first we create a set of sequence regions
	region_pair_vector segments;
	for (fragment_ptr_list::iterator j = ext.segments.begin(); j != ext.segments.end(); ++j)
		segments.push_back( make_region_pair( target, query, **j ) );

	if( use_viterbi_training ) aligner.constrained_viterbi_train( segments.begin(), segments.end() );
	else trainer.constrained_bw_train( segments.begin(), segments.end(), num_threads );
}

template<class EXTENDER, class ALIGNER, class TRAINER>
bool least<EXTENDER,ALIGNER,TRAINER>::align( dna_sequence_ptr query, extension &ext ) {
	using namespace boost;
	using namespace std;

	ostream &out = *output_stream;

	// first we create a set of sequence regions
	region_pair_vector segments;
	for (fragment_ptr_list::iterator i = ext.segments.begin(); i != ext.segments.end(); ++i)
		segments.push_back( make_region_pair( target, query, **i ) );


	annotation_ptr ann;
	dna_alignment_sequence_ptr alignmenta(new dna_alignment_sequence()), alignmentb(new dna_alignment_sequence());
	pairwise_dna_alignment final(alignmenta, alignmentb, 0.0);

	aligner.set_interface_size(40);
	tie(final, ann) = aligner.constrained_align_and_annotate(segments.begin(),segments.end());

	// Check viterbi cut threshold.
	if( use_viterbi_cut && final.score < viterbi_cut_threshold ) return false;

	// Don't output alignments if we are training.
	if( training_mode ) return true;

	// Output MAF format.
	out << "a score=" << final.score << endl;

	out << "s " << target->tags["accession"] << '\t' << (*ext.segments.front())[0].a;
	int len = ((*ext.segments.back())[0].b - (*ext.segments.front())[0].a + 1);
	out << ' ' << len;
	out << ' ' << target->tags["strand"];
	out << ' ' << target->data.size() << ' ';
	for( dna_alignment_sequence_data::const_iterator i = final.a->data.begin(); i != final.a->data.end(); ++i ) out << *i;
	out << endl;

	int real_len = 0;
	for( dna_alignment_sequence_data::const_iterator i = final.a->data.begin(); i != final.a->data.end(); ++i ) if( *i != dna_alignment_alpha::GAP) real_len++;



	out << "s " << query->tags["accession"] << '\t' << (*ext.segments.front())[1].a;
	out << ' ' << ((*ext.segments.back())[1].b - (*ext.segments.front())[1].a + 1);
	out << ' ' << query->tags["strand"];
	out << ' ' << query->data.size() << ' ';
	for( dna_alignment_sequence_data::const_iterator i = final.b->data.begin(); i != final.b->data.end(); ++i ) out << *i;
	out << endl;

	// Output extra least information.
	out << "#feast_model";
	for( annotation::const_iterator i = ann->begin(); i != ann->end(); ++i ) out << ' ' << i->t << ' ' << i->size();
	out << endl;

	out << "#feast_extension_score " << ext.score.as_base() << endl << endl;

	if( real_len != len ) {
		cerr << "LENGTH MISMATCH real " << real_len << " claimed " << len << endl;
		exit(1);
	}
	return true;
}

template<class EXTENDER, class ALIGNER, class TRAINER>
boost::shared_ptr< std::vector<point> > least<EXTENDER,ALIGNER,TRAINER>::find_hit_points(dna_sequence_ptr query) {
	using namespace std;
	using namespace __gnu_cxx;
	using namespace boost;
	vector<dna_sequence_region_ptr> sequences;
	dna_sequence_region_ptr seqA( new dna_sequence_region( target ) );
	sequences.push_back( seqA );
	dna_sequence_region_ptr seqB( new dna_sequence_region( query ) );
	sequences.push_back( seqB );

	/*fragment allf;
	allf.add(closed_interval(0,in.a.seq->data.size()-1));
	allf.add(closed_interval(0,in.b.seq->data.size()-1));
	seed_parameters p2 = computeMSP(in,allf);
	cerr << "Minimum seed weight " << get<0>(p2) << endl;*/

	/*
	seed_parameters p1 = computeOSP(in, f, 2000000);
	seed_parameters p2 = computeMSP(in,f);
	seed_parameters p =  get<0>(p1) > get<0>(p2) ? p1 : p2;
	random_seed_factory sf(get<1> (p));
	seed myseed = sf.create(get<0> (p));
	cerr << "seed " << myseed.pattern() << " " << get<0>(p) << " " << get<1>(p) << " " << get<2>(p2) << endl;
*/

	vector<seed> seeds;

	// non-coding seeds from "Designing multiple simultaneous seeds for DNA similarity search"
	seeds.push_back(seed("111100110011111"));
	seeds.push_back(seed("1101110000001100001111"));
	seeds.push_back(seed("1111000110000000110111"));

	// coding seeds from "Designing multiple simultaneous seeds for DNA similarity search"
	seeds.push_back(seed("1101100001101101101"));
	seeds.push_back(seed("1101101101100000001101"));
	seeds.push_back(seed("1111000011000011001011"));

	typedef fast_hash_hit_data<dna_alpha> hit_data;

	shared_ptr< vector<point> > hit_points( new vector<point>() );

	point_hit_builder< back_insert_iterator< vector<point> > > point_builder( back_inserter(*hit_points) );

	for( vector<seed>::iterator i = seeds.begin(); i != seeds.end(); ++i ) {
#ifdef DEBUG_OUTPUT
		cerr << "seed " << i->pattern_string() << ' ' << i->weight() << ' ' << i->length() << endl;
#endif
		fast_hash_factory< hit_data >hasher(*i, true);
		boost::shared_ptr< hit_data > sequence_hash = hasher.create_from_ptr(sequences.begin(),sequences.end());
		sequence_hash->enumerate( point_builder );
	}
	return hit_points;
}


/*


transformed_fragment_ptr_vector_ptr least<EXTENDER,ALIGNER,TRAINER>::remove_overlaps( transformed_fragment_ptr_vector &frags) {
	using namespace std;
	transformed_fragment_ptr_vector_ptr newfrags( new transformed_fragment_ptr_vector());
	sort(frags.begin(), frags.end(), fragment_algorithm::tfptr_yxlt());
	for (transformed_fragment_ptr_vector::iterator i = frags.begin(), next =	++frags.begin();;) {
		if (next != frags.end() && (**next).a().y == (**i).a().y	&& (**next).a().x <= (**i).b().x)
			++next;
		newfrags->push_back(*i);

		if (next == frags.end()) break;

		i = next;
		++next;

	}
	return newfrags;
}

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_min_extend_score( double s ) { min_extend_score = s; }

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_min_chain_score( double s ) { min_chain_score.set_base(s); }

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_anchor_search_step( int s ) { ext_step = s; }

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_anchor_search_max_box( int s ) {
	ext_max = s;
	hcost = -1.0/(double)ext_max;
}
template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_training_mode( bool t ) { training_mode = t; }


template<class EXTENDER, class ALIGNER, class TRAINER>
least<EXTENDER,ALIGNER,TRAINER>::seed_parameters least<EXTENDER,ALIGNER,TRAINER>::computeOSP( input_data &in, fragment &f, int random_sample ) {
	using namespace std;

	// estimate region composition
	blitz::TinyVector<double, 4> q = composition_frequency(*(in.a.stats), f[0]);
	q += composition_frequency(*(in.b.stats), f[1]);
	q /= (double) (f[0].size() + f[1].size());

	// compute match probability for a column
	double pr_match = dot(q, q);
	double area = (double) f[0].size() * (double) f[1].size();
	unsigned int w = (unsigned int) (log((double) random_sample / area) / log(pr_match));
	unsigned int l = (unsigned int) ((double) w * 1.4);
	if (l < 5)	l = 5; // minimum length for random seeds is 5
	if (w < 3)	w = 3; // minimum weight 3
	unsigned int expected_random_hits = (unsigned int) (pow(pr_match,	(double) w) * area);
	return seed_parameters(w, l, expected_random_hits);
}

template<class EXTENDER, class ALIGNER, class TRAINER>
least<EXTENDER,ALIGNER,TRAINER>::seed_parameters least<EXTENDER,ALIGNER,TRAINER>::computeMSP( input_data &in, fragment &f ) {
	using namespace std;

	// estimate region composition
	blitz::TinyVector<double, 4> q = composition_frequency(*(in.a.stats), f[0]);
	q += composition_frequency(*(in.b.stats), f[1]);
	q /= (double) (f[0].size() + f[1].size());

	// compute match probability for a column
	double pr_match = dot(q, q);
	double a = -1.0/hcost, b = -1.0/(hcost+vcost);
	double area = b*b+(a-b)*b;

	unsigned int w = (unsigned int) (log((double) 0.5 / area) / log(pr_match));
	unsigned int l = (unsigned int) ((double) w * 1.4);
	if (l < 5)	l = 5; // minimum length for random seeds is 5
	if (w < 3)	w = 3; // minimum weight 3
	unsigned int expected_random_hits = (unsigned int) (pow(pr_match,	(double) w) * area);
	return seed_parameters(w, l, expected_random_hits);
}

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_hcost( double s ) {	hcost = s; }

template<class EXTENDER, class ALIGNER, class TRAINER>
void least<EXTENDER,ALIGNER,TRAINER>::set_vcost( double s ) {	vcost = s; }
*/

#endif
