/** This just passes the integer through. */

#ifndef NULL_HASH_H_
#define NULL_HASH_H_

template < typename HashKey >
struct null_hash {
   size_t operator()( HashKey const &x ) const { return x; }
};

#endif /* DJB2_HASH_H_ */
