#ifndef __SEEDFACTORY__
#define __SEEDFACTORY__

#include <string>
#include <boost/shared_ptr.hpp>
#include "seed.h"

class seed_factory {
	protected:
		virtual seed make_seed( unsigned int weight ) = 0;

	public:
		seed_factory() {}
		virtual ~seed_factory() {}

      inline seed create( unsigned int weight ) { 
         return make_seed( weight );
      }  
};

typedef boost::shared_ptr<seed_factory> seed_factory_ptr;

#endif
