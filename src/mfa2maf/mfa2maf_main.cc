#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <btl/fasta_reader.h>
#include <btl/mdna_iupac_gap.h>
#include <btl/sequence.h>

#include "utility.h"

// forward declare version variable
extern const char* build_string;

int main( int argc, char* argv[] ) {
	   // load program options
   using namespace boost::program_options;
   using namespace std;

   string program_name("mfa2maf build ");
   program_name += string(build_string);

    options_description desc("Allowed options");
   desc.add_options()
      ( "help", "produce help message" )
      ( "in,i", value<string>(), "input mfa file" )
      ( "out,o", value<string>(), "output maf file" )
   ;
   variables_map vm;
   store( parse_command_line( argc, argv, desc ), vm );
   notify(vm);

   if( (argc < 2) || vm.count("help")) {
      cout << program_name << endl << desc << endl;
      return 1;
   }

   require_option( vm, "in", 2 );

   // open mfa file
   ifstream in_file;
   in_file.open( vm["in"].as<string>().c_str(), ios::in );
   if( in_file.fail() ) {
      cerr << "unable to open " << vm["in"].as<string>() << " for reading." << endl;
      return 3;
   }

   // default to cout
   ostream *out = &cout;

   // open output if specified
   ofstream out_file;
   if( vm.count("out") ) {
      out_file.open( vm["out"].as<string>().c_str(), ios::out|ios::trunc );
      if( out_file.fail() ) {
         cerr << "unable to open " << vm["out"].as<string>() << " for writing." << endl;
         return 4;
      }
      out = &out_file;
   }

   // Output MAF format, we assume a single alignment.
	*out << "##maf version=1" << endl;
   *out << "a" << endl;

   // input sequence type
   typedef btl::sequence<btl::mdna_iupac_gap> input_sequence;
   input_sequence inseq;

   while(1) {
      in_file >> btl::fasta_reader( inseq );
		if( in_file.fail() ) break;

      // count gaps
      int gaps = count( inseq.data.begin(), inseq.data.end(), btl::mdna_iupac_gap::GAP );
      int original_size = inseq.data.size() - gaps;
 
      *out << "s " << inseq.tags["accession"] << "\t0\t" << original_size << "\t+\t" << original_size << "\t";
		copy( inseq.data.begin(), inseq.data.end(), ostream_iterator<btl::mdna_iupac_gap_symbol>(*out) );
      *out << endl;
   }

   in_file.close();
   if( vm.count("out") ) out_file.close();
   
   return EXIT_SUCCESS;

}