#ifndef __BTL_BASIC_TREE__
#define __BTL_BASIC_TREE__

#include <string>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

//--------------------------------------------------------------------------
// This file defines a basic phylogenetic tree using an adjacency_list.
// For most purposes this graph definition should be fine.
//--------------------------------------------------------------------------
namespace btl {
   typedef boost::adjacency_list< 
      boost::vecS, 
      boost::vecS, 
      boost::bidirectionalS,
      boost::property<boost::vertex_name_t, std::string>, 
      boost::property<boost::edge_weight_t, double> > basic_tree;

   typedef boost::property_map<basic_tree, boost::vertex_index_t>::type basic_index_map;
   typedef boost::property_map<basic_tree, boost::vertex_name_t>::type basic_vertex_name;
   typedef boost::property_map<basic_tree, boost::edge_weight_t>::type basic_edge_weight;
   
   typedef boost::graph_traits<basic_tree> basic_traits;
   typedef basic_traits::vertex_descriptor basic_vertex;
   typedef basic_traits::vertex_iterator basic_vertex_iter;
   typedef basic_traits::edge_iterator basic_edge_iter;
   typedef basic_traits::in_edge_iterator basic_in_edge_iter;
   typedef basic_traits::out_edge_iterator basic_out_edge_iter;
}

#endif
