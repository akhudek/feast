#ifndef __BTL_ANNOTATED_SEQUENCE__
#define __BTL_ANNOTATED_SEQUENCE__

#include <vector>
#include <set>
#include "btl/sequence.h"
#include "btl/annotation.h"

namespace btl {
	template<
		typename ALPHABET, 
		typename CONTAINER = stl::vector<ALPHABET>,
		typename ANNOTATION_CONTAINER = stl::set<btl::annotation> >
	class annotated_sequence : public sequence<ALPHABET,CONTAINER> {
		public:
			ANNOTATION_CONTAINER = annotations;
			
		private:
			inline void _copy( annotated_sequence<ALPHABET,CONTAINER,ANNOTATION_CONTAINER> const &seq ) {
				annotations = seq.annotations;
				sequence<ALPHABET,CONTAINER>::_copy( seq );
			}
		public:
			/// Default constructor.
			annotated_sequence() {}

         /// Copy constructor.
         annotated_sequence( annotated_sequence<ALPHABET,CONTAINER,ANNOTATION_CONTAINER> const &seq ) {
            this->_copy( seq );
         }

         /// Clear all data.
         void clear() {
            *this = sequence<ALPHABET,CONTAINER>();
         }

	};
}
#endif /*__BTL_ANNOTATED_SEQUENCE__*/
