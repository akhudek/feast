#ifndef __BTL_PHYLOGENETICTREE__
#define __BTL_PHYLOGENETICTREE__

#include <btl/basic_tree.h>

namespace btl {
   template<class Graph = basic_tree>
   class phylogenetic_tree {
      public:
         typedef boost::graph_traits<Graph> traits;
         typedef typename traits::vertex_descriptor vertex;
         typedef Graph graph;

         vertex root;
         graph g;
   };
}
#endif
