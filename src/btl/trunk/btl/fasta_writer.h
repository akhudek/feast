#ifndef BTL_FASTA_WRITER__
#define BTL_FASTA_WRITER__

#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <boost/spirit/include/classic_core.hpp>
#include <btl/sequence.h>

namespace btl {
   class fasta_writer {
	private:
		// base non-templated class containing templated class interface
		struct base {
			virtual void write( std::ostream &output ) const = 0;
			virtual base *clone() const = 0;
			virtual ~base() {}
		};
		
		//-------------------------------------------------------------------
		// templated inner class, all public because no one should see this
		//-------------------------------------------------------------------
		template<typename ALPHABET,  typename CONTAINER >
		struct fasta_writer_t : public base {
			sequence<ALPHABET, CONTAINER> const &seq;
			
			// templated class constructor
			fasta_writer_t( sequence<ALPHABET, CONTAINER> const &newSeq ) : seq( newSeq ) {}
			
			// clone this class
			base *clone() const {
				return new fasta_writer_t<ALPHABET, CONTAINER>(seq);
			}
			
			// Read seqeunce from stream into sequence object
			void write( std::ostream &output ) const {
				using namespace std;
				typename sequence<ALPHABET,CONTAINER>::tag_map_const_iter tag;
				
				// write header
				output << '>';
				
				// write accession
				tag = seq.tags.find("accession");
				if( tag == seq.tags.end() ) {
					output << "sequence";
				} else { 
					output << tag->second;
				}
				
				// write description
				tag = seq.tags.find("description");
				if( tag != seq.tags.end() ) {
					output << ' ' << tag->second;
				}
				
				// finish header
				output << endl;
				
				// write sequence, 60 characters per line
				int cnt = 0;
				for( typename CONTAINER::const_iterator i = seq.data.begin(); i != seq.data.end(); ++i ) {
					output << *i;
					
					// track characters per line
					cnt++;
					if( cnt == 60 ) {
						output << endl;
						cnt = 0;
					}
				}
				// finish last line if we need to
				if( cnt != 0 ) {
					output << endl;
				}
			}
		};
		
		// pointer to templated function
		base *bp;
		
	public:
		// non-templated class constructor, uses type deduction to avoid
		// gunk when calling 
		template<typename ALPHABET,typename CONTAINER >
		fasta_writer( sequence<ALPHABET,CONTAINER> const &seq ) {
			bp = new fasta_writer_t<ALPHABET,CONTAINER>( seq );
		}
		
		// copy constructor
		fasta_writer( fasta_writer const &other ) {
			bp = other.bp->clone();
		}
		
		// copy operator
		void operator=( fasta_writer const &other ) {
			delete bp;
			bp = other.bp;
		}
		
		// destructor
		~fasta_writer() {
			delete bp;
		}
		
		friend std::ostream &operator<<( std::ostream &out, fasta_writer const &fobj ) {
			fobj.bp->write( out );
			return out;
		}
   }; // end class fasta_format
}

#endif
