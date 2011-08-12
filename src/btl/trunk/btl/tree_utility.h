#ifndef __BTL_TREE_UTILITY__
#define __BTL_TREE_UTILITY__

#include <utility>
#include <vector>
#include <cassert>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

namespace btl {
   // get the parent vertex discriptor
   template<typename Vertex, typename Graph>
   inline Vertex parent( Vertex v, Graph &g ) {
      using namespace boost;
      return source( parent_edge( v, g ), g );
   }

   template<typename Vertex, class Graph>
   inline typename boost::graph_traits<Graph>::edge_descriptor 
      parent_edge( Vertex v, Graph &g ) {
      using namespace boost;
      typedef typename graph_traits<Graph>::in_edge_iterator in_edge_iter;
      
      int indegree = in_degree(v, g);
      assert( indegree == 1 );         // ensure we are not root

      std::pair< in_edge_iter, in_edge_iter > iep = in_edges(v, g);
      return *iep.first;
   }


   template<typename Vertex, typename Graph>
   inline Vertex add_child( Vertex v, Graph &g ) {
   	  using namespace boost;
      Vertex child = add_vertex(g);
      add_edge( v, child, g );
      return child;
   }

   template<typename Vertex, class Graph>
   inline std::vector<Vertex> children( Vertex v, Graph &g ) {
      using namespace boost;
      typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iter;
      std::vector<Vertex> results;
      results.reserve( out_degree(v, g) );

      std::pair<out_edge_iter, out_edge_iter> oep;
      for( oep = out_edges(v, g); oep.first != oep.second; ++oep.first ) {
         results.push_back(target(*oep.first, g));
      }
      return results;
   }

   template<typename Vertex, class Graph>
   inline bool is_leaf( Vertex v, Graph &g ) {
      using namespace boost;
      return out_degree(v, g) == 0;
   }

   template<typename Vertex, class PhylogeneticTree>
   inline bool is_root( Vertex v, PhylogeneticTree &t ) {
      return v == t.root;
   }
};

#endif
