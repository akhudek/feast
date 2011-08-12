#include "entropy.h"

#include <iostream> 
#include <stdexcept>
#include <ctime>

// Constructor.
entropy::entropy() {
   using namespace std;
#ifndef WINDOWS
   devRandom.open("/dev/urandom", ios::in | ios::binary );
   if( devRandom.fail() ) throw runtime_error( "entroy() : cannot open /dev/urandom" );
#endif
}

// Deconstructor.
entropy::~entropy() {
#ifndef WINDOWS
   //devRandom.close();
#endif
}

// Return an integer.
int entropy::get_int() {
#ifndef WINDOWS
   int val = devRandom.get();
   if( devRandom.fail() || devRandom.bad() ) throw std::runtime_error( "entropy::get_int() : error reading from /dev/urandom" );
   return val;

#else
   // No /dev/urandom on Windows, use time.
   return (int) time(NULL);
#endif
}

