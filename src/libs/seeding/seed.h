/**
 * Provides storage and a few methods to get seed properties.
 */

#ifndef SEED_H__
#define SEED_H__
#include <string>
#include <vector>
#include <algorithm>

class seed {
public:
	/** Data type for seed mask. */
	typedef std::vector<bool> seed_mask;

protected:
	/** Seed pattern. */
	std::string p;

	/** Seed pattern as mask. */
	seed_mask m;

public:
	/** Default constructor. Not intended to be used. */
	seed();

	/** Seed constructor. */
	seed(std::string const &pattern);

	/** Return the length of a seed. */
	unsigned int length() const;

	/** computes the self overlap */
	int self_overlap() const;

	/** Returns a const reference to the pattern. */
	std::string const &pattern_string() const;

	/** Returns a const reference to the pattern mask. */
	seed_mask const &pattern_mask() const;

	/** Returns the number of 1's in the seed.
	 *
	 *  The current implementation is O(L), where L is the length of the seed.
	 */
	unsigned int weight() const;
};

// Inline function implementations.

inline std::string const &seed::pattern_string() const {
	return p;
}

inline seed::seed_mask const &seed::pattern_mask() const {
	return m;
}

inline unsigned int seed::weight() const {
	return std::count(m.begin(), m.end(), true );
}

inline unsigned int seed::length() const {
	return p.size();
}

#endif
