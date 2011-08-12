#ifndef __BTL_NEWICK_WRITER__
#define __BTL_NEWICK_WRITER__

#include <iostream>
#include <string>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <btl/tree_utility.h>

namespace btl {
   class newick_writer {
      private:
         // base non-templated class containing templated class interface
         struct base {
            virtual void write( std::ostream &output ) const = 0;
            virtual base *clone() const = 0;
            virtual ~base() {}
         };

         //-------------------------------------------------------------------
         // templated inner class, all private because no one should see this
         //-------------------------------------------------------------------
         template < class PTree, typename NTrait = boost::vertex_name_t, typename DTrait = boost::edge_weight_t>
         struct newick_writer_t : public base {
            PTree const &tree;

            // templated class constructor
            newick_writer_t( PTree const &t ) : tree( t ) {}

            // clone this class
            base *clone() const {
               return new newick_writer_t<PTree,NTrait,DTrait>(tree);
            }
 
            // Read seqeunce from stream into sequence object
            void write( std::ostream &output ) const {
               using namespace std;
               newick_output<NTrait,DTrait> newick_out(output);
               boost::depth_first_search( tree.g, visitor(newick_out).root_vertex(tree.root) );
            }
         };

         // pointer to templated function
         base *bp;

         // newick_output class: writes tree via depth first search
         template<typename NameTrait, typename DistanceTrait>
         class newick_output : public boost::default_dfs_visitor {
            private:
               std::ostream &output;

            public:
               newick_output(std::ostream &out) : output(out) {}

               template <typename Vertex, class Graph>
               void discover_vertex( Vertex v, Graph &g ) {
                  using namespace std;
                  using namespace boost;
                  typedef typename property_map<Graph, vertex_index_t>::type index_map;
                  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iter;
                  typedef pair<out_edge_iter,out_edge_iter> out_edge_iter_pair;
                  
                  index_map index = get(vertex_index, g);

                  // output ',' if we need to
                  if( index[v] != 0 ) {
                     Vertex p = parent(v, g);
                     out_edge_iter_pair eip = out_edges(p,g);
                     if( target(*eip.first,g) != v ) {
                        output << ", ";
                     }
                  }
                  
                  if( !is_leaf(v, g) ) {
                     output << '(';
                  }
              }

               template <typename Vertex, class Graph>
               void finish_vertex( Vertex v, Graph &g ) { 
                  using namespace std;
                  using namespace boost;
                  typedef typename property_map<Graph, vertex_index_t>::type index_map;
                  typedef typename graph_traits<Graph>::edge_descriptor edge;

                  NameTrait name_trait;
                  DistanceTrait distance_trait;
                  index_map index = get(vertex_index, g);
                     
                  if( !is_leaf(v, g) ) {
                     output << ')';
                  }

                  string name = get( name_trait, g, v );
                  if( name.size() > 0 ) {
                    output << name;
                  }

                  // if we are not root 
                  if( index[v] != 0 ) {
                     edge e = parent_edge(v, g);
                     double dist = get( distance_trait, g, e );
                     output << ':' << dist;

                  } else {
                     output << ';';
                  }
               }
         };

      public:
         // non-templated class constructor, uses type deduction to avoid gunk when calling 
         template<typename NTrait, class PTree>
         newick_writer( PTree const &g ) {
            bp = new newick_writer_t<PTree,NTrait>(g);
         }
         
         template<class PTree>
         newick_writer( PTree const &g ) {
            bp = new newick_writer_t<PTree>(g);
         }
         
         // copy constructor
         newick_writer( newick_writer const &other ) {
            bp = other.bp->clone();
         }

         // copy operator
         void operator=( newick_writer const &other ) {
            delete bp;
            bp = other.bp;
         }

         // destructor
         ~newick_writer() {
            delete bp;
         }
    
         friend std::ostream &operator<<( 
            std::ostream &out, 
            newick_writer const &fobj ) {
            fobj.bp->write( out );
            return out;
         }
   }; // end class newick_format
}

#endif
