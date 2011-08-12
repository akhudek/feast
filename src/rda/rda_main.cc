#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <btl/fasta_writer.h>

#include "ereal.h"
#include "dna_sequence.h"
#include "main_io.h"
#include "io.h"

#include "dna_alignment_sequence.h"
#include "rda_functions.h"
#include "krda_improve.h"
#include "needleman.h"
#include "nw_model_parameters.h"
#include "alignment_functions.h"
#include "utility.h"
#include "annotation.h"

#define PARAMETER_ERROR 1

// forward declare version variable
extern const char* build_string;

int main( int argc, char* argv[] ) {
	   // load program options
   using namespace boost::program_options;
   using namespace std;
   
   string program_name("rda build ");
   program_name += string(build_string);

   options_description desc("Allowed options");
   desc.add_options()
      ( "help", "produce help message" )
      ( "target,t", value<string>(), "first sequence" )
      ( "query,q", value<string>(), "second sequence" )
      ( "out,o", value<string>(), "output alignment file" )
      ( "parameter-file,p", value<string>(), "set parameters from file" )
      ( "align-method", value<unsigned int>()->default_value(0), "alignment method: 0 RDA, 1 needleman" )
      ( "block-size,k", value<unsigned int>()->default_value(40), "maximum block size" )
      ( "block-breaks", "add pure gap sites as block breaks" )
      ( "output-maf", "output alignment in MAF format")
   ;
   variables_map vm;
   store( parse_command_line( argc, argv, desc ), vm );
   notify(vm);

   if( (argc < 2) || vm.count("help")) {
      cout << program_name << endl << desc << endl;
      return 1;
   }

   require_option( vm, "target", PARAMETER_ERROR );
   require_option( vm, "query", PARAMETER_ERROR );

   unsigned int K = vm.count("block-size") ? vm["block-size"].as<unsigned int>() : 40;
     
   // precompute lookup tables
   ereal::init();

   try {
      dna_sequence_ptr target,query;
		target = load_sequence( vm["target"].as<string>() );
      query = load_sequence( vm["query"].as<string>() );

      require_size_above( *target, 1 );
      require_size_above( *query, 1 );

      krda_improve<needleman> krda;
      krda.set_k(K);
      needleman nw;
		
      // load parameter file if needed
      if( vm.count("parameter-file") ) {
         nw_model_parameters_ptr mp = load_parameters_from_file( vm["parameter-file"].as<string>() );
			krda.set_parameters( *mp );
         nw.set_s( mp->pr_open, mp->p, mp->q );
         nw.set_mean_gap_length( mp->mean_gap_length );
      }

      if( vm.count("block-breaks") )  krda.set_block_marker( true );

      // build final alignment by aligning each identified region with needleman wunch
      dna_alignment_sequence_ptr alignmenta = new_dna_alignment_sequence();
      dna_alignment_sequence_ptr alignmentb = new_dna_alignment_sequence();
      pairwise_dna_alignment final = pairwise_dna_alignment( alignmenta, alignmentb, 0 );
      dna_sequence_region_ptr alla = new_dna_sequence_region(target), allb = new_dna_sequence_region(query);
		
		vector< pair<dna_sequence_region_ptr,dna_sequence_region_ptr> > training_set;
		training_set.push_back( make_pair(alla,allb) );

      // open output file now so that if there is a problem we don't waste time aligning stuff
      ofstream out_file;
      if( vm.count("out") ) {
         out_file.open( vm["out"].as<string>().c_str(), ios::out|ios::trunc );
         if( out_file.fail() ) {
            cerr << "unable to open " << vm["out"].as<string>() << " for writing" << endl;
            return 4;
         }
      }

		annotation_ptr myann;
      switch( vm["align-method"].as<unsigned int>() ) { 
         case 1:
            final = nw.align(*alla,*allb);
            break;
         case 0:
            final = krda.align(*alla,*allb);
            break;
         default:
            throw runtime_error("unknown alignment method");
      }
      
      // label alignments
      final.a->tags["accession"] = target->tags["accession"];
      final.b->tags["accession"] = query->tags["accession"];
      string info = program_name;
      switch( vm["align-method"].as<unsigned int>() ) { 
         case 1:
            info += string(" nm score=") + boost::lexical_cast<string>((double)final.score);
            break;
         case 0:
            info += string(" rda k=") + boost::lexical_cast<string>(K);
            info += string(" score=") + boost::lexical_cast<string>((double)final.score);
      }
      final.a->tags["description"] = final.b->tags["description"] = info;
     
      ostream *out = vm.count("out") ? &out_file :  &cout;
		
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
   } catch( exception &e ) {
      cerr << "FATAL: " << e.what() << endl;
   }

   return 0;
}
