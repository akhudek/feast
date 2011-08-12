#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <btl/fasta_writer.h>
#include <btl/logspace.h>

#include "dna_sequence.h"
#include "main_io.h"
#include "io.h"
#include "cape_pairwise_aligner.h"
#include "nw_model_parameters.h"
#include "annotation_io.h"
#include "utility.h"

#define SUCCESS 0
#define PARAMETER_ERROR 1
#define RUNTIME_ERROR 2

// forward declare version variable
extern const char* build_string;

int main( int argc, char* argv[] ) {
   // load program options
   using namespace boost::program_options;
   using namespace std;
   
   string program_name("cape build ");
   program_name += string(build_string);

#ifdef VITERBI_EXTENSIONS
	program_name += " VITERBI_EXTENSIONS";
#endif

   options_description required_parameters("Required parameters");
   required_parameters.add_options()
      ( "target,t", value<string>(), "first sequence" )
      ( "query,q", value<string>(), "second sequence" )
   ;

   options_description optional_parameters("Optional parameters") ;
   optional_parameters.add_options()
      ( "help", "produce help message" )
      ( "out,o", value<string>(), "output alignment file" )
//      ( "train", "train parameters with Baum-Welch, output is a set of trained parameters" )
      ( "write-annotation,a", value<string>(), "write region annotation to file in GFF format" )
      ( "parameters,p", value<string>(), "set parameters from file" )
      ( "initialize,i", value<string>(), "set parameters from initialization file")
      ( "threads", value< int >()->default_value(1), "number of threads to use" )
      ( "extension-limit", value< int >()->default_value(200), "set anchor extension limit" )
      ( "random-sample", value< int >()->default_value(500), "set number of random anchors to aim for" )
      ( "required-region-size", value< unsigned int >()->default_value(2000), "set the maximum region size between anchors" )
//		( "inference-points", value< unsigned int >()->default_value(4), "set the number of inference points to use. Set to zero to use supplied parameters." )
//		( "max-training-passes", value< unsigned int >()->default_value(1000), "the maximum number of training passes to use")
		( "chain-filter-distance", value< unsigned int >()->default_value(0), "distance from longest chain within which to keep fragments. Set to 0 to disable.")
		( "greedy-threshold", value< double >()->default_value(70.0), "set the threshold under which we switch to greedy anchoring" )
		( "enable-random-model", "enable the random model" )
		( "output-maf", "output alignment in MAF format")
   ;

   options_description cmd_line_parameters;
   cmd_line_parameters.add(required_parameters);
   cmd_line_parameters.add(optional_parameters);

   variables_map vm;
   store( parse_command_line( argc, argv, cmd_line_parameters ), vm );
   notify(vm);

   if( (argc < 2) || vm.count("help")) {
      cout << program_name << endl << cmd_line_parameters << endl;
      return SUCCESS;
   }
   require_option( vm, "target", PARAMETER_ERROR );
   require_option( vm, "query", PARAMETER_ERROR );
   require_one_option( vm, "parameters", "initialize", PARAMETER_ERROR );

   if( vm["threads"].as<int>() < 1 ) {
      cerr << "Must have at least one thread.";
		return PARAMETER_ERROR;
	}

	// precompute lookup tables
	ereal::init();

   try {
      // Check for number of models and choose specialized model if we are using one model.
   	int nmodels = 0;

      if( vm.count("parameters") ) {
         ifstream infile(vm["parameters"].as<string> ().c_str());
         if (infile.fail()) throw runtime_error("set_parameters_from_file: could not open parameter file");
         // extract number of models, skip random parameters
         nmodels = extract_number_of_models(infile);
      } else if( vm.count("initialize") )  {
         ifstream infile(vm["initialize"].as<string> ().c_str());
         if (infile.fail()) throw runtime_error("initialize_parameters_from_file: could not open initialization file");
         // extract models from initialization format
         nmodels = read<int>(infile);
      }

      dna_sequence_ptr target,query;
		target = load_sequence( vm["target"].as<string>() );
      query = load_sequence( vm["query"].as<string>() );

      require_size_above( *target, 1 );
      require_size_above( *query, 1 );

      // Measure background frequencies of input sequences.
      blitz::TinyVector<unsigned int,4> countT,count;
      sequence_statistics_ptr target_stats = measure_statistics( *target );
      countT = composition_frequency(*target_stats,closed_interval(0,target->data.size()-1));
      sequence_statistics_ptr query_stats = measure_statistics( *query );
      count = composition_frequency(*target_stats,closed_interval(0,target->data.size()-1));
      count += countT;
      unsigned int total = blitz::sum(count);
      blitz::TinyVector<double,4> q;
      for(int i = 0; i < 4; ++i ) q(i) = (double)count(i)/(double)(total);
      
		cape_pairwise_aligner cape( nmodels );
      cape.set_background(q);
      cape.set_extension_limit( vm["extension-limit"].as<int>() );
      cape.set_random_sample( (double) vm["random-sample"].as<int>() );
      cape.set_required_region_size( vm["required-region-size"].as<unsigned int>() );
//		if( vm["inference-points"].as<unsigned int>() > 0 ) cape.set_active_parameter_sets( vm["inference-points"].as<unsigned int>() );
		cape.set_custom_scoring(true);
		cape.set_threshold( vm["greedy-threshold"].as<double>() );
      cape.set_num_threads( vm["threads"].as<int>() );
		//cape.set_training_passes( vm["max-training-passes"].as<unsigned int>() );
		cape.set_chain_filter_distance( (int) vm["chain-filter-distance"].as<unsigned int>() );
		if( !vm.count("enable-random-model") ) cape.mma_disable_random_model();
		if( vm.count("parameters") ) cape.set_parameters_from_file(vm["parameters"].as<string>());
		else cape.initialize_parameters_from_file(vm["initialize"].as<string>());
      
      // open output file now so that if there is a problem we don't waste time aligning stuff
      ofstream out_file;
      if( vm.count("out") ) {
         out_file.open( vm["out"].as<string>().c_str(), ios::out|ios::trunc );
         if( out_file.fail() ) {
            cerr << "unable to open " << vm["out"].as<string>() << " for writing" << endl;
            return 4;
         }
      }

		// open annotation output file now so that if there is a problem we don't waste time aligning stuff
      ofstream aout_file;
      if( vm.count("write-annotation") ) {
         aout_file.open( vm["write-annotation"].as<string>().c_str(), ios::out|ios::trunc );
         if( aout_file.fail() ) {
            cerr << "Unable to open " << vm["write-annotation"].as<string>() << " for writing." << endl;
            throw runtime_error("Could not open file.");
         }
      }

		if( vm.count("train") ) {
			if( vm.count("out") ) {
				cape.train( target, query ).write_parameters( out_file );
				out_file.close();
			} else {
				cape.train( target, query ).write_parameters( cout );
			}

		} else {
			// do alignment
			pairwise_dna_alignment final;
			annotation_ptr ann;
			boost::tie( final, ann ) = cape.align_and_annotate( target, query );

			// label alignments
			final.a->tags["accession"] = target->tags["accession"];
			final.b->tags["accession"] = query->tags["accession"];
			string info = program_name;
			info += string(" score=") + boost::lexical_cast<string>( final.score );
			final.a->tags["description"] = final.b->tags["description"] = info;

			ostream *out = &cout;

			if( vm.count("out") ) out = &out_file;

			if( vm.count("output-maf") ) {
				// Output MAF format.
				*out << "##maf version=1" << endl;
				*out << "a score=" << final.score << endl;
				*out << "s " << final.a->tags["accession"] << "\t0\t" << target->data.size() << "\t+\t" << target->data.size() << "\t";
				for( dna_alignment_sequence_data::const_iterator j = final.a->data.begin(); j != final.a->data.end(); ++j ) *out << *j;
				*out << endl;
				*out << "s " << final.b->tags["accession"] << "\t0\t" << query->data.size() << "\t+\t" << query->data.size() << "\t";
				for( dna_alignment_sequence_data::const_iterator j = final.b->data.begin(); j != final.b->data.end(); ++j ) *out << *j;
				*out << endl;
			} else {
				*out << btl::fasta_writer( *final.a );
				*out << btl::fasta_writer( *final.b );
			}

			if( vm.count("out") ) out_file.close();

			if( vm.count("write-annotation") ) write_seqa_gff_annotation( aout_file, *ann, final );
		}
   } catch( exception &e ) {
      cerr << "FATAL: " << e.what() << endl;
   }

   return 0;
}
