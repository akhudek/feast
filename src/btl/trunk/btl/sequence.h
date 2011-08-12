#ifndef __BTL_SEQUENCE_
#define __BTL_SEQUENCE_

#include <vector>
#include <string>
#include <map>

namespace btl {
   	template<typename ALPHABET, typename CONTAINER = std::vector<typename ALPHABET::symbol> >
   	class sequence {
	public:
		typedef ALPHABET alphabet;
       	typedef CONTAINER data_type;

         /// The seqeunce data is stored here. Data is encoded by ALPHABET.
         data_type data;

         /// Property map to hold sequence information.
         typedef std::map<std::string,std::string> tag_map;
         typedef tag_map::iterator tag_map_iter;
         typedef tag_map::const_iterator tag_map_const_iter;
         tag_map tags;

		private:
			inline void _copy( sequence<ALPHABET,CONTAINER> const &seq ) {
            data = seq.data;
            tags = seq.tags;
			}

		public:
         /// Default constructor.
         sequence() {}

         /// Copy constructor.
         sequence( sequence<ALPHABET,CONTAINER> const &seq ) {
         		this->_copy(seq);
         }

         /// Clear all data.
         void clear() {
            *this = sequence<ALPHABET,CONTAINER>();
         }
   };

}

#endif
