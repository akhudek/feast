#ifndef BTL_FASTA_READER__
#define BTL_FASTA_READER__

#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <boost/spirit/include/classic_core.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <btl/sequence.h>

namespace btl {
   class fasta_reader {
	private:
		// base non-templated class containing templated class interface
		struct base {
			virtual void read( std::istream &input ) const = 0;
			virtual base *clone() const = 0;
			virtual ~base() {}
		};

		//-------------------------------------------------------------------
		// templated inner class, all public because no one should see this
		//-------------------------------------------------------------------
		template< typename ALPHABET, typename CONTAINER >
		struct fasta_reader_t : public base {
			sequence<ALPHABET, CONTAINER> &seq;

			// templated class constructor
			fasta_reader_t( sequence<ALPHABET, CONTAINER> &newSeq ) : seq( newSeq ) {}

			// clone this class
			base *clone() const {
				return new fasta_reader_t<ALPHABET, CONTAINER>(seq);
			}

			// Read sequence from stream into sequence object
			void read( std::istream &input ) const {
				using namespace std;
				using namespace boost::spirit::classic;
				string instr;

				// clear sequence data
				seq.data.clear();
				seq.tags.clear();

				// find next sequence
				input >> skipws;
				getline( input, instr, '\n' );
				boost::algorithm::trim(instr);

				// parse header (set accession and description
				string accession, description;
				bool fully_parsed = parse( instr.begin(), instr.end(),
						// Fasta header grammar.
						(
							'>'
							>> (+(print_p - space_p))[assign_a(accession)]
							>> !( +space_p >> +print_p)[assign_a(description)]
						)
					).full;

				if( !fully_parsed ) {
					// clear tags and set failed bit
					seq.tags.clear();
					input.setstate( ios::failbit );
					return;
				}
				boost::algorithm::trim(description);
				seq.tags["accession"] = accession;
				seq.tags["description"] = description;

				// read sequence
				while( 1 ) {
					// skip white space
					input >> ws;
					if( input.eof() || ( input.peek() == '>' )) break;

					// read rest of sequence line skipping whitespace
					getline( input, instr, '\n' );
					stringstream parseline(stringstream::in|stringstream::out);
					parseline << instr;
					while( 1 ) {
						parseline >> instr;
						if( parseline.fail() ) break;
						copy( instr.begin(), instr.end(), back_inserter(seq.data) );
					}
				}
			}
		};

		// pointer to templated function
		base *bp;

	public:
		// non-templated class constructor, uses type deduction to avoid
		// gunk when calling
		template<typename ALPHABET, typename CONTAINER  >
		fasta_reader( sequence<ALPHABET,CONTAINER> &seq ) {
			bp = new fasta_reader_t<ALPHABET,CONTAINER>( seq );
		}

		// copy constructor
		fasta_reader( fasta_reader const &other ) {
			bp = other.bp->clone();
		}

		// copy operator
		void operator=( fasta_reader const &other ) {
			delete bp;
			bp = other.bp;
		}

		// destructor
		~fasta_reader() {
			delete bp;
		}

		friend std::istream &operator>>( std::istream &in, fasta_reader const &fobj ) {
			fobj.bp->read( in );

			return in;
		}
   }; // end class fasta_format
}

#endif
