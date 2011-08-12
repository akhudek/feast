#ifndef __BTL_ANNOTATION__
#define __BTL_ANNOTATION__
// BTL annotation class
// 
// An annotation is a range associated with a sequence.


namespace btl {
	class annotation {
	public:
		enum annotation_type { repeat = 0 };
		
		int start;
		int end;
		annotation_type type;
		
		// default constructor
		annotation() : start(0), end(0), type(repeat) {}
	};
}

#endif /*__BTL_ANNOTATION__*/
