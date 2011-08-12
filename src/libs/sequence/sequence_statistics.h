#ifndef SEQUENCE_STATISTICS_H__
#define SEQUENCE_STATISTICS_H__
#include <vector>
#include <boost/shared_ptr.hpp>
#include "triplet_table.h"

// for now, we only have one entry, eventually we may have more 
struct sequence_statistics {
   // number of repeat bases up to point i
   std::vector<unsigned int> repeat_bases;   

   // composition up to point i
   std::vector<unsigned int> A;
   std::vector<unsigned int> T;
   std::vector<unsigned int> C;
   std::vector<unsigned int> G;

   // measure triplets
   std::vector< triplet_table<unsigned int> > tri;
};

typedef boost::shared_ptr<sequence_statistics> sequence_statistics_ptr;


#endif
