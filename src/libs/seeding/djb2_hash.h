/** djb2 hash function.
 * from http://www.cse.yorku.ca/~oz/hash.html
 * this algorithm (k=33) was first reported by dan bernstein many years ago in comp.lang.c. another version of
 * this algorithm (now favored by bernstein) uses xor: hash(i) = hash(i - 1) * 33 ^ str[i]; the magic of number
 * 33 (why it works better than many other constants, prime or not) has never been adequately explained.
 */

#ifndef DJB2_HASH_H_
#define DJB2_HASH_H_

template < typename Container >
struct djb2_hash {
   size_t operator()( Container const &x ) const {
      unsigned long hash = 5381;
      hash = 5381;
      for( typename Container::const_iterator i = x.begin(); i != x.end(); ++i ) {
		hash = ((hash << 5) + hash) ^ (int)*i;
      }
      return hash;
   }
};

#endif /* DJB2_HASH_H_ */
