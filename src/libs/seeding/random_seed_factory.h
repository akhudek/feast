#ifndef __RANDOMSEEDFACTORY__
#define __RANDOMSEEDFACTORY__

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/shared_ptr.hpp>
#include <string>
#include <vector>

#include "seed_factory.h"

#define RSF_MAX_L 41

class random_seed_factory : public seed_factory {
	private:
      typedef std::vector<unsigned int> overlap_vector;
		typedef boost::mt19937 base_generator;
		typedef boost::shared_ptr<base_generator> base_generator_ptr;
		typedef boost::uniform_01<base_generator> uniform_generator;
		typedef boost::shared_ptr<uniform_generator> uniform_generator_ptr;

		// static array of all data
		static unsigned int const overlap_values[RSF_MAX_L][RSF_MAX_L];
		
		// overlap cutoff
      overlap_vector overlap_cutoff;
      unsigned int length;

		// random number generators
      base_generator_ptr baseg;
      uniform_generator_ptr uniformg;

      // generate a random seed string
      std::string generate_seed_string( unsigned int w );

   protected:
      seed make_seed( unsigned int weight );

   public:
      random_seed_factory( unsigned int length );
};

// simple allocator
seed_factory_ptr new_random_seed_factory( unsigned int L );

#endif
