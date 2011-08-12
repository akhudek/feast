/*
 *  explore_main.cc
 *  rda
 *
 *  Created by Alexander K. Hudek on 2008-09-16.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <btl/gff_feature.h>
#include <btl/algorithm.h>
#include <cmath>
#include <fstream>
#include <functional>
#include <algorithm>
#include "explore_io.h"
#include "match_data.h"

//#define DEBUG
// forward declare version variable
extern const char* build_string;

int main(int argc, char* argv[]) {
	// load program options
	using namespace boost::program_options;
	using namespace std;

	string program_name("explore build ");
	program_name += string(build_string);

	string feast_str("feast");
	string global_str("global");

	options_description required_parameters("Required parameters");
	required_parameters.add_options()
	("alignment,a", value<std::string>(), "alignment")
	("reference,r", value<std::string>(), "reference sequence")
	("features,f", value<std::string>(), "features (relative to reference)")
	("input-type,t", value<std::string>()->default_value(global_str), "alignment type")
	;

	options_description optional_parameters("Optional parameters");
	optional_parameters.add_options()
	("help", "produce help message")("window",value<unsigned int>()->default_value(50), "window length")
	("conservation-threshold", value<double>()->default_value(0.75), "level to consider feature conserved")
	("conservation-min-length", value<int>()->default_value(50), "minimum length of positions >= threshold to consider feature conserved")
	("output-conserved", value<std::string>(),"file to append conserved features to")
	("output-match-data", value<std::string>(), "write match data to file")
	("alignment-coverage", "measure alignment coverage")
	;

	options_description cmd_line_parameters;
	cmd_line_parameters.add(required_parameters);
	cmd_line_parameters.add(optional_parameters);

	variables_map vm;
	store(parse_command_line(argc, argv, cmd_line_parameters), vm);
	notify(vm);

	if ((argc < 2) || vm.count("help")) {
		cout << program_name << endl << cmd_line_parameters << endl;
		return 1;
	}
	if (!vm.count("alignment")) {
		cout << "missing alignment" << endl;
		return 2;
	}

	if (!vm.count("reference")) {
		cout << "missing reference sequence" << endl;
		return 2;
	}

	if (vm["input-type"].as<std::string> () != feast_str
			&& vm["input-type"].as<std::string> () != global_str) {
		cout << "input-type must be either 'global' or 'feast'" << endl;
		return 4;
	}

	int window = (int) vm["window"].as<unsigned int> ();

	try {
		// Load reference sequence.
		ex_sequence_ptr reference = load_sequence(vm["reference"].as<std::string>());

		// Load alignment.
		ex_sequence_ptr_vector_ptr alignment = load_alignment(vm["alignment"].as<std::string>());

		if (alignment->size() < 2) {
			cerr << "Alignment has less than two sequences." << endl;
			exit(3);
		}

		if (alignment->size() > 2 && vm["input-type"].as<std::string>() != feast_str) {
			cerr << "Global alignment type does not support more than two sequences." << endl;
			exit(3);
		}

		// create match data
		match_data_ptr md(new match_data(reference->data.size(), 0.0));

		unsigned int next_i = 0;
		ex_sequence_ptr A = (*alignment)[next_i++], B = (*alignment)[next_i++], RA, OA;
		while (1) {
			// Find reference sequence.
			string idA = A->tags["accession"], idB = B->tags["accession"];

			// Remove pair id from feast output.
			if (vm["input-type"].as<std::string> () == feast_str) {
				string::size_type i = idA.rfind('_'), j = idB.rfind('_');
				if (i == string::npos || j == string::npos) {
					cerr << "Pair not feast output: " << endl;
					cerr << idA << endl;
					cerr << idB << endl;
					exit(5);
				}
				idA.erase(i);
				idB.erase(j);
			}

			if (reference->tags["accession"] != idA) {
				if (reference->tags["accession"] != idB) {
					cerr << "Cannot identify reference in pair: " << endl;
					cerr << idA << endl;
					cerr << idB << endl;
					exit(5);
				}
				RA = B;
				OA = A;
			} else {
				RA = A;
				OA = B;
			}

			// Default to entire reference range.
			closed_interval r(0, reference->data.size() - 1);

			// Adjust for feast.
			if (vm["input-type"].as<std::string> () == feast_str) {
				stringstream descRef(RA->tags["description"]);
				descRef >> r.a;
				descRef >> r.b;
			}

			// filter alignment according to reference sequence
			ex_sequence_ptr_vector_ptr falignment(new ex_sequence_ptr_vector);
			ex_sequence_ptr FA(new ex_sequence()), FB(new ex_sequence());
			falignment->push_back(FA);
			falignment->push_back(FB);
			for (unsigned int c = 0; c < RA->data.size(); ++c) {
				if (RA->data[c] != btl::mdna_iupac_gap::GAP) {
					FA->data.push_back(RA->data[c]);
					FB->data.push_back(OA->data[c]);
				}
			}

			if( vm.count("alignment-coverage") ) fill_match_data(*md, falignment, r, window, coverage() );
			else fill_match_data(*md, falignment, r, window, max_function()  );

			// move to next pair
			if (next_i == alignment->size())
				break;
			assert( alignment->size() - next_i >= 2 );
			A = (*alignment)[next_i++], B = (*alignment)[next_i++];
		}

		if (vm.count("output-match-data")) {
			fstream match_file;
			string match_filename = vm["output-match-data"].as<std::string> ();
			match_file.open(match_filename.c_str(), fstream::app | fstream::out);

			for (int i = 0; i < (int) md->size(); ++i)
				match_file << (*md)[i] << endl;

			match_file.close();
		}

		if( vm.count("alignment-coverage") ) exit(0);

		// count frequencies
		int resolution = 20;
		vector<int> rfreq(resolution, 0);
		for (int i = 0; i < (int) md->size(); ++i) {
			rfreq[floor((*md)[i] * (resolution - 1))]++;
		}

		if (vm.count("features")) {
			// load features
			vector<btl::gff_feature> features;
			string feat_filename = vm["features"].as<std::string> ();
			fstream feature_file(feat_filename.c_str(), fstream::in);
			btl::read_gff_feature_set(feature_file, back_inserter(features));
			feature_file.close();

			fstream conserved_file;

			if (vm.count("output-conserved")) {
				string conserved_filename = vm["output-conserved"].as<std::string> ();
				conserved_file.open(conserved_filename.c_str(), fstream::app|fstream::out);
			}

			// tally features conserved
			double T = vm["conservation-threshold"].as<double> ();
			int L = vm["conservation-min-length"].as<int> ();

			map<string, int> conserved_features;
			map<string, int> total_features;
			map<string, bool> feature_set;
			for (vector<btl::gff_feature>::const_iterator i = features.begin(); i != features.end(); ++i) {
				string feature_key;
				feature_key = i->feature + "-"
						+ vm["reference"].as<std::string> () + "-"
						+ boost::lexical_cast<string>(i->start) + "-"
						+ boost::lexical_cast<string>(i->end);

				// Skip feature if we already explored it. This happens due to alternative spliced genes.
				if (feature_set.find(feature_key) != feature_set.end())
					continue;
				else
					feature_set[feature_key] = true;

				if (total_features.find(i->feature) == total_features.end()) {
					total_features[i->feature] = 0;
					conserved_features[i->feature] = 0;
				}
				total_features[i->feature] += 1;

				// check for conservation
				if (i->start < (int) md->size() && (i->end - window + 1) < (int)md->size()) {
					int ctally = 0;
					for (int j = i->start; j <= (i->end - window + 1); ++j) {
						if ((*md)[j] >= T) ctally += 1;
					}
					if (ctally >= L) {
						conserved_features[i->feature] += 1;
						if (conserved_file.good())
							conserved_file << i->feature << "-" << vm["reference"].as<std::string> ()
									<< "-" << i->start << "-" << i->end << endl;
					}
				}
			}

			// output features conserved
			for (map<string, int>::const_iterator i = total_features.begin(); i != total_features.end(); ++i)
				cout << i->first << "\t" << conserved_features[i->first]<< "\t" << i->second << endl;
		}

		// output total positions
		cout << reference->data.size() << endl;

		// output frequencies
		for (int i = 0; i < 20; ++i) {
			cout << rfreq[i] << endl;
		}

	} catch (exception &e) {
		cerr << "FATAL: " << e.what() << endl;
	}

	return 0;
}
