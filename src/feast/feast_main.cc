#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdio>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <btl/fasta_writer.h>
#include <btl/logspace.h>

#include "dna_sequence.h"
#include "io.h"
#include "main_io.h"
#include "multi_alignment_model.h"
#include "mm_alignment.h"
#include "sm_extension.h"
#include "least.h"
#include "utility.h"

#define PARAMETER_ERROR 1

// forward declare version variable
extern const char* build_string;

std::string program_name("feast build ");

const char * masked_warning = "WARNING: %s is all lowercase and will be considered completely masked.\n"
                        "         Use --ignore-mask to change this behaviour.\n";

// Set parameters from vm.
template <class LALIGNER>
void set_parameters( LALIGNER &laligner, boost::program_options::variables_map &vm ) {
	using namespace std;

	// We initialize the target first so that initialized parameters can
	// use the background frequency.
	dna_sequence_ptr target = load_sequence(vm["target"].as<std::string> ());
   if( vm.count("ignore-mask") ) unmask( *target );
   else if( is_all_masked( *target ) ) fprintf( stderr, masked_warning, vm["target"].as<std::string>().c_str() );
	target->tags["strand"] = "+";

   require_size_above( *target, 1 );
	laligner.set_target( target );

	if( vm.count("parameters") ) laligner.set_parameters_from_file(vm["parameters"].as<std::string>());
	else if( vm.count("initialize") ) laligner.initialize_parameters_from_file(vm["initialize"].as<std::string>());

	laligner.set_drop_threshold( - vm["drop-threshold"].as<double>() );
	laligner.set_ungapped_drop_threshold( - vm["ungapped-drop-threshold"].as<double>() );
	laligner.set_cut_threshold( vm["cut-threshold"].as<double>() );
	laligner.set_ungapped_cut_threshold( vm["ungapped-cut-threshold"].as<double>() );
	laligner.set_num_threads( vm["threads"].as<unsigned int>() );
	laligner.set_stop_at_repeat( vm.count("end-at-repeat") );

	if( vm.count("viterbi-cut-threshold") )	laligner.set_viterbi_cut( vm["viterbi-cut-threshold"].as<double>() );
	if( vm.count("use-viterbi-training") ) laligner.set_viterbi_training( true );

}

// Templated class to work with templated least class.
template <class LALIGNER>
void do_alignment( LALIGNER &laligner, boost::program_options::variables_map &vm ) {
	using namespace std;

	set_parameters(laligner,vm);

	ofstream seq_out;
	ostream *out = &cout;
	if( vm.count("out") ) {
		seq_out.open(vm["out"].as<std::string>().c_str(), ios::trunc|ios::out );
		if( seq_out.fail() ) throw runtime_error("could not open sequence output file");
		out = &seq_out;
	}
	laligner.set_output_stream(out);

	dna_sequence_ptr_vector_ptr query_set = load_sequence_set(vm["query"].as<string>() );
  bool lowercase_query = true;
   for( dna_sequence_ptr_vector::iterator i = query_set->begin(); i != query_set->end(); ++i ) {
      if( vm.count("ignore-mask") ) unmask( **i );
      else if( !is_all_masked( **i ) ) lowercase_query = false;
   }
   if( lowercase_query ) fprintf( stderr, masked_warning, vm["query"].as<std::string>().c_str() );
   
	// Start MAF output.
	*out << "##maf version=1 program=feast" << endl;
	*out << "# -t " << vm["target"].as<string>() << " -q " << vm["query"].as<string>();
	if( vm.count("parameters") ) *out << " -p " << vm["parameters"].as<string>();
	else *out << " -i " << vm["initialize"].as<string>();
	if( vm.count("positive-strand") ) *out << " --positive-strand";
	if( vm.count("end-at-repeat") ) *out << " --end-at-repeat";
  	if( vm.count("ignore-mask") ) *out << " --ignore-mask";
	*out << " --drop-threshold " << vm["drop-threshold"].as<double>();
	*out << " --ungapped-drop-threshold " << vm["ungapped-drop-threshold"].as<double>();
	*out << " --cut-threshold " << vm["cut-threshold"].as<double>();
	*out << " --ungapped-cut-threshold " << vm["ungapped-cut-threshold"].as<double>();
	*out << endl;
	*out << "# " << program_name << endl;

	for( dna_sequence_ptr_vector::iterator i = query_set->begin(); i != query_set->end(); ++i ) {
		(**i).tags["strand"] = "+";
		laligner.search( *i );

		if( !vm.count("positive-strand") ) {
			// Build reverse compliment of query.
			dna_sequence_ptr query_rev(new dna_sequence());
			query_rev->tags = (**i).tags;
			for( dna_sequence::data_type::reverse_iterator j = (**i).data.rbegin(); j != (**i).data.rend(); ++j )
				query_rev->data.push_back( j->compliment() );

			query_rev->tags["strand"] = "-";
			laligner.search(query_rev);
		}
	}
}

#define DIVIDER "---------------------------------"

// Templated class to work with templated least class.
template <class LALIGNER, class TRAINER>
void do_training( LALIGNER &laligner, boost::program_options::variables_map &vm, TRAINER &trainer ) {
	using namespace std;

	set_parameters(laligner,vm);

	ofstream out;
	if( !vm.count("out") ) {
		cerr << "Standard output unsupported for training." << endl;
		exit(1);
	}

	dna_sequence_ptr_vector_ptr query_set = load_sequence_set(vm["query"].as<string>() );
   bool lowercase_query = true;
   for( dna_sequence_ptr_vector::iterator i = query_set->begin(); i != query_set->end(); ++i ) {
      if( vm.count("ignore-mask") ) unmask( **i );
      else if( !is_all_masked( **i ) ) lowercase_query = false;
   }
   if( lowercase_query ) fprintf( stderr, masked_warning, vm["query"].as<std::string>().c_str() );

   cerr << "INITIAL PARAMETERS" << endl;
	trainer.pretty_print_parameters(cerr);
	cerr << DIVIDER << endl << endl;

	double stop_delta = vm["stop-delta"].as<double>();
	unsigned int pass = 1;
	ereal last_pr = 0.0;
	while(1) {

		cerr << "Pass " << pass << endl << endl;

		laligner.start_training();
		for( dna_sequence_ptr_vector::iterator i = query_set->begin(); i != query_set->end(); ++i ) {
			(**i).tags["strand"] = "+";
			laligner.search(*i);

			if( !vm.count("positive-strand") ) {
				// Build reverse compliment of query.
				dna_sequence_ptr query_rev(new dna_sequence());
				query_rev->tags = (**i).tags;
				for( dna_sequence::data_type::reverse_iterator j = (**i).data.rbegin(); j != (**i).data.rend(); ++j )
					query_rev->data.push_back( j->compliment() );

				query_rev->tags["strand"] = "-";
				laligner.search(query_rev);
			}
		}

		ereal pr = laligner.stop_training();


		trainer.pretty_print_parameters( cerr );
		cerr << "Pr:    " << pr << endl;
		double delta = last_pr.as_base() - pr.as_base();
		cerr << "Delta: " << delta << endl;

		out.open(vm["out"].as<std::string>().c_str(), ios::trunc|ios::out );
		if( out.fail() ) throw runtime_error("could not open sequence output file");
		trainer.write_parameters(out);
		out.close();

		cerr << DIVIDER << endl << endl;

		++pass;
		last_pr = pr;

		if( abs(delta) < stop_delta ) break;
	}
	cerr << "TRAINING COMPLETE." << endl;
}


int main(int argc, char* argv[]) {
	// load program options
	using namespace boost::program_options;
	using namespace std;


	program_name += string(build_string);

#ifdef VITERBI_EXTENSIONS
	program_name += " VITERBI_EXTENSIONS";
#endif

#ifdef ANCHOR_AT_ORIGIN
	program_name += " ANCHOR_AT_ORIGIN";
#endif

	options_description required_parameters("Required parameters");
	required_parameters.add_options()
		("target,t", value<std::string> (),"single target sequence in fasta format")
		("query,q", value<std::string> (),"multiple query sequences in fasta format")
		;

	options_description model_parameters("Model parameters (one required)");
	model_parameters.add_options()
		("parameters,p", value<std::string> (), "set parameters from file")
		("initialize,i", value<std::string>(), "set parameters from initialization file")
		;

	options_description extension_parameters("Extension parameters (optional)");
	extension_parameters.add_options()
		("drop-threshold", value<double>()->default_value(65.0), "drop threshold for extension algorithm")
		("cut-threshold", value<double>()->default_value(20.0), "cut off threshold for extension algorithm")
		("viterbi-cut-threshold", value<double>(), "filter out alignments based on viterbi score")
		("ungapped-drop-threshold", value<double>()->default_value(10.0), "drop threshold for ungapped extension algorithm")
		("ungapped-cut-threshold", value<double>()->default_value(15.0), "cut off threshold for ungapped extension algorithm")
		("positive-strand", "only explore positive strand")
		("end-at-repeat", "stop extension when a repeat is encountered")
      ("ignore-mask", "treat input files as all uppercase (no masking")
		;

	options_description training_parameters("Training parameters (optional)");
	training_parameters.add_options()
		("train", "train parameters with Baum-Welch, output is a set of trained parameters" )
		("stop-delta", value<double>()->default_value(5.0), "probability change threshold at which to stop training" )
		("use-viterbi-training", "use Viterbi training instead of Baum-Welch")
		;

	options_description optional_parameters("Miscellaneous parameters (optional)");
	optional_parameters.add_options()
		("help", "produce help message")
		("out,o",value<std::string> (), "redirect standard output to file")
		("threads", value<unsigned int>()->default_value(1), "maximum number of threads")
		;


	options_description cmd_line_parameters;
	cmd_line_parameters.add(required_parameters);
	cmd_line_parameters.add(model_parameters);
	cmd_line_parameters.add(extension_parameters);
	cmd_line_parameters.add(training_parameters);
	cmd_line_parameters.add(optional_parameters);

	variables_map vm;
	store(parse_command_line(argc, argv, cmd_line_parameters), vm);
	notify(vm);

	if ((argc < 3) || vm.count("help")) {
		cout << program_name << endl << cmd_line_parameters << endl;
		return 1;
	}
	require_option(vm,"target",PARAMETER_ERROR);
	require_option(vm,"query",PARAMETER_ERROR);
	if( vm["threads"].as<unsigned int>() == 0 ) {
		cout << "Must have at least one thread." << endl;
		exit(1);
	}

	//cerr << program_name << endl;

	// precompute lookup tables
	ereal::init();

	try {
		// Check for number of models and choose specialized model if we are using one model.
		int nmodels = 0;

		if( vm.count("parameters") && vm.count("initialize") ) {
			cerr << "Must specify either -p or -i, but not both." << endl;
			exit(1);
		} else if( vm.count("parameters") ) {
			ifstream infile(vm["parameters"].as<std::string> ().c_str());
			if (infile.fail()) throw runtime_error("set_parameters_from_file: could not open parameter file");
			// extract number of models, skip random parameters
			nmodels = extract_number_of_models(infile);
		} else if( vm.count("initialize") )  {
			ifstream infile(vm["initialize"].as<std::string> ().c_str());
			if (infile.fail()) throw runtime_error("initialize_parameters_from_file: could not open initialization file");
			// extract models from initialization format
			nmodels = read<int>(infile);
		} else {
			cerr << "One of -p or -i is required." << endl;
			exit(1);
		}

		sm_extension single_model;
		multi_alignment_model multi_model(nmodels);
		mm_alignment aligner(nmodels);

		// Always use multi_alignment_model for training so that trained parameters are used
		// by the extender after each pass.
		if( nmodels == 1 && !vm.count("train") ) {
			least<sm_extension,mm_alignment,multi_alignment_model> myleast(single_model,aligner,multi_model);
			do_alignment(myleast, vm);
		} else {
			least<multi_alignment_model,mm_alignment,multi_alignment_model> myleast(multi_model,aligner,multi_model);
			if( !vm.count("train") ) do_alignment(myleast, vm);
			else do_training(myleast,vm,multi_model);
		}

	} catch (exception &e) {
		cerr << "FATAL: " << e.what() << endl;
	}

	return 0;
}
