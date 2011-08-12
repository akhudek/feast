#ifndef SIMILARITY_REGION_H__
#define SIMILARITY_REGION_H__
#include <deque>
#include <utility>
#include <boost/shared_ptr.hpp>
#include "closed_interval.h"

typedef std::pair<closed_interval,closed_interval> similarity_region;

typedef std::deque<similarity_region> similarity_region_deque;
typedef similarity_region_deque::iterator similarity_region_deque_iter;
typedef similarity_region_deque::const_iterator similarity_region_deque_citer;

typedef boost::shared_ptr<similarity_region_deque> similarity_region_deque_ptr;

// allocator helper fucntion
inline similarity_region_deque_ptr new_similarity_region_deque() {
   similarity_region_deque_ptr result( new similarity_region_deque() );
   return result;
}
#endif
