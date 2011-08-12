#ifndef __BTL_BINARY_TREE_UTILITY__
#define __BTL_BINARY_TREE_UTILITY__

#include <utility>
#include <vector>
#include <cassert>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

namespace btl {
   template<typename Vertex, class Graph>
   inline Vertex left_child( Vertex v, Graph &g ) {
      using namespace boost;
      typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iter;

      // we need at least one edge
      //if( out_degree(v, g) < 1 ) return graph_traits<GRAPH>::null_vertex();
     if( out_degree(v, g) < 1 ) return -1;

      std::pair< out_edge_iter, out_edge_iter > oep = out_edges(v,g);
      return target(*oep.first,g);
   }

   template<typename Vertex, class Graph>
   inline Vertex right_child( Vertex v, Graph &g ) {
      using namespace boost;
      typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iter;

      // we need at least two edges
      //if( out_degree(v, g) < 2 ) return graph_traits<GRAPH>::null_vertex();
     if( out_degree(v, g) < 2 ) return -1;

      std::pair< out_edge_iter, out_edge_iter > oep = out_edges(v,g);
      return target(*(++oep.first),g);
   }

}

#endif
