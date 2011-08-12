#ifndef CHAR_HASH_FUNCTION_H_
#define CHAR_HASH_FUNCTION_H_

// Hash function for generic characters.

// From sdbm:
//
// polynomial conversion ignoring overflows
// hash[i] = hash[i-1] * 65599 + c[i]
// 
// 65587 is even better supposedly
template < typename Container >
struct hashf {
   size_t operator()( Container const & x ) const {
      unsigned long hash = 0;
      unsigned int c = 0;

      for( typename Container::const_iterator i = x.begin(); i != x.end(); i++ ) {
         // take one character
         c = *i;
         hash = c + (hash << 6) + (hash << 16) - hash;
      }
      return hash;
   }
};

#endif /*CHAR_HASH_FUNCTION_H_*/
