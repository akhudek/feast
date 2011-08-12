#include "fragment.h"

std::ostream &operator<<( std::ostream &out, fragment const &f ) {
   out << "[";
   for( std::vector<closed_interval>::const_iterator i = f.regions.begin(); i != f.regions.end(); i++ ) {
      out << " " << *i;
   }
   out << " ]";
   return out;
}

bool operator<<( fragment const &fa, fragment const &fb ) {
   typedef std::vector<closed_interval>::const_iterator region_citer;
   for( region_citer i = fa.regions.begin(), j = fb.regions.begin(); i != fa.regions.end(); i++, j++ ) {
      if( (*i).a <= (*j).b ) return false;
   }
   return true;
}

bool operator>>( fragment const &fa, fragment const &fb ) {
   typedef std::vector<closed_interval>::const_iterator region_citer;
   for( region_citer i = fa.regions.begin(), j = fb.regions.begin(); i != fa.regions.end(); i++, j++ ) {
      if( (*i).b >= (*j).a ) return false;
   }
   return true;
}
