#ifndef ENTROPY_H__
#define ENTROPY_H__
#define WINDOWS

// This really should be a singleton.
#include <fstream>

class entropy {
   protected:
      std::fstream devRandom;

   public:
      entropy();
      ~entropy();

      int get_int();
};

#endif
