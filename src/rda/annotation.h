#ifndef ANNOTATION_H__
#define ANNOTATION_H__
#include "annotation_interval.h"
#include <boost/shared_ptr.hpp>
#include <list>

// annotation
typedef std::list<annotation_interval> annotation;
typedef annotation::iterator annotation_iter;
typedef annotation::reverse_iterator annotation_riter;
typedef annotation::const_iterator annotation_citer;
typedef boost::shared_ptr<annotation> annotation_ptr;

std::ostream &operator<<( std::ostream &out, annotation const &a );

#endif

