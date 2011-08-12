#ifndef __BTL_NEWICK_READER__
#define __BTL_NEWICK_READER__

#include <string>
#include <cstdlib>
#include <btl/tree_utility.h>
#define BOOST_SPIRIT_RULE_SCANNERTYPE_LIMIT 2
#include <boost/spirit/classic_core.hpp>
#include <boost/spirit/tree/parse_tree.hpp>
#include <boost/spirit/tree/ast.hpp>

namespace btl {
   class newick_reader {
      private:
         // grammar class for newick format
         struct newick_grammar : public boost::spirit::grammar<newick_grammar> {
            enum nodeID { INNER_NODE = 1, LEAF_NODE, LIST, QUOTED_LABEL, UNQUOTED_LABEL, QUOTED_CHAR, UNQUOTED_CHAR, DISTANCE };
                                                                                
            template <typename ScannerT>
            struct definition {
               typedef char const* iterator_t;
               typedef typename ScannerT::iteration_policy_t iteration_policy_t;
               typedef typename ScannerT::action_policy_t action_policy_t;
              
               // define pt_nodes
               typedef boost::spirit::pt_match_policy<iterator_t> pt_match_policy_t;
               typedef boost::spirit::scanner_policies <iteration_policy_t, pt_match_policy_t> pt_scanner_policies_t;
               typedef boost::spirit::scanner <iterator_t, pt_scanner_policies_t> pt_scanner_t;
               
               // define ast_nodes
               typedef boost::spirit::ast_match_policy<iterator_t> ast_match_policy_t;
               typedef boost::spirit::scanner_policies <iteration_policy_t, ast_match_policy_t> ast_scanner_policies_t;
               typedef boost::spirit::scanner <iterator_t, ast_scanner_policies_t> ast_scanner_t;
               
               // define no_nodes
               typedef boost::spirit::match_policy match_policy_t;
               typedef boost::spirit::scanner_policies <iteration_policy_t, match_policy_t, action_policy_t> policies_t;
               typedef boost::spirit::scanner <iterator_t, policies_t> non_tree_scanner_t;
               typedef boost::spirit::parser_context<> parser_context_t;
         
               // define rules
               boost::spirit::rule<ast_scanner_t> label, subtree;
         
               // set some id's for special nodes
               boost::spirit::rule< ast_scanner_t, parser_context_t, boost::spirit::parser_tag<QUOTED_LABEL> > quoted_label;
               boost::spirit::rule< ast_scanner_t, parser_context_t, boost::spirit::parser_tag<QUOTED_CHAR> > quoted_print_p;
               boost::spirit::rule< ast_scanner_t, parser_context_t, boost::spirit::parser_tag<UNQUOTED_LABEL> > unquoted_label;
               boost::spirit::rule< ast_scanner_t, parser_context_t, boost::spirit::parser_tag<UNQUOTED_CHAR> > unquoted_print_p;
               boost::spirit::rule< ast_scanner_t, parser_context_t, boost::spirit::parser_tag<DISTANCE> > branch_length;
               boost::spirit::rule< pt_scanner_t, parser_context_t, boost::spirit::parser_tag<LIST> > descendant_list;
               boost::spirit::rule< pt_scanner_t, parser_context_t, boost::spirit::parser_tag<LEAF_NODE> > leaf_node;
               boost::spirit::rule< pt_scanner_t, parser_context_t, boost::spirit::parser_tag<INNER_NODE> > inner_node;
               boost::spirit::rule<non_tree_scanner_t> whitespace;
                                                                                
               boost::spirit::rule<ScannerT> tree;
                                                                                
               definition( newick_grammar const& self ) {
                  using namespace boost::spirit;

                  // NOTE: Stay away from infix_node_d and inner_node_d. Use no_node_d instead. I've had no end to problems
                  // with the former.
                                                                                
                  // define some characters character
                  chlit<char> quote('\'');
                  chlit<char> lbrac('(');
                  chlit<char> rbrac(')');
                  chlit<char> ldelim(',');
                  chlit<char> ndelim(':');
                  chlit<char> tree_end(';');
                                                                                
                  whitespace = *(ch_p(' ') | ch_p('\t') | ch_p('\n') | ch_p('\r') ) ;
                                                                                
                  // define a quoted label
                  quoted_print_p = (print_p - quote) | ( no_node_d[ quote ] >> quote );
                  quoted_label = inner_node_d[ quote >> leaf_node_d[ +quoted_print_p ] >> quote ];
                                                                                
                  // define an unquoted label
                  unquoted_print_p = print_p - quote - space_p - lbrac
                     - rbrac - ch_p('[') - ch_p(']') - ndelim - tree_end
                     - ldelim;
                  unquoted_label = leaf_node_d[ +unquoted_print_p ];
                                                                                
                  // define a regular label
                  label = no_node_d[ whitespace ] >> ( quoted_label | unquoted_label ) >> no_node_d[ whitespace ];
                                                                                
                  // define a branch length
                  branch_length =  no_node_d[ whitespace ]
                                >> leaf_node_d[ !ch_p('-') >> +digit_p
                          >> !( ch_p('.') >> +digit_p ) 
                          >> !( ( ch_p('e') | ch_p('E') ) 
                                >> ( ch_p('-') | ch_p('+') ) >> +digit_p ) ]
                          >> no_node_d[ whitespace ];
                                                                                
                  // define a leaf
                  leaf_node = gen_ast_node_d[ label >> !( no_node_d[ ndelim ] >> branch_length ) ];
                  
                  // define a basic list
                  descendant_list =  no_node_d[ whitespace ] 
                           >> gen_ast_node_d[ no_node_d[ lbrac ] >> subtree >> *( no_node_d[ ldelim ] >> subtree ) >> no_node_d[ rbrac ] ] 
                           >> no_node_d[ whitespace ];
                                                                                
                  // inner node, can have optional label and distance
                  inner_node = descendant_list >> gen_ast_node_d[ !label >> !( no_node_d[ ndelim ] >> branch_length ) ];
                                                                                
                  // a subtree is either another list or a leaf
                  subtree = gen_pt_node_d[ inner_node | leaf_node ];
            
                  // the full tree
                  tree = gen_pt_node_d[ inner_node ] >> no_node_d[ tree_end ];
                  
               }
                                                                                
               // start of grammar
               boost::spirit::rule<ScannerT> const& start() const { 
                  return tree; 
               }
            };
         };

         // base non-templated class containing templated class interface
         struct base {
            virtual void read( std::istream &input ) const = 0;
            virtual base *clone() const = 0;
            virtual ~base() {}
         };

         //-------------------------------------------------------------------
         // templated inner class, all public because no one should see this
         //-------------------------------------------------------------------
         template< typename PTree >
         struct newick_reader_t : public base {
            PTree &tree;
            typedef boost::graph_traits<typename PTree::graph> graph_traits;
            typedef typename graph_traits::vertex_descriptor vertex_t;

            typedef boost::spirit::tree_match<char const *> tree_match_t;
            typedef typename tree_match_t::tree_iterator tree_iterator_t;
            typedef typename tree_match_t::container_t tree_container_t;
            typedef typename tree_match_t::node_t tree_node_t;

            void set_values( vertex_t nodeId, tree_iterator_t field, tree_iterator_t end ) const {
               using namespace boost;
               using namespace std;
               // check for LABEL
               if( (field->value.id() == newick_grammar::QUOTED_LABEL)
                   || (field->value.id() == newick_grammar::UNQUOTED_LABEL)
                   || (field->value.id() == newick_grammar::QUOTED_CHAR)
                   || (field->value.id() == newick_grammar::UNQUOTED_CHAR) ) {

                  string name( field->value.begin(), field->value.end() );
                  put( vertex_name, tree.g, nodeId, name );
                  field++;
               } 

               // no value attached so return
              if( field == end ) return;
                                                                                
               // check for distance
               if( field->value.id() == newick_grammar::DISTANCE ) {
                  string textValue( field->value.begin(), field->value.end() );
                  double numericValue = atof( textValue.c_str() );
                  typename graph_traits::edge_descriptor pedge;
                  pedge = parent_edge( nodeId, tree.g );
                  put( edge_weight, tree.g, pedge, numericValue );
               }
            }

            void assign_tree( vertex_t nodeId, tree_node_t &node ) const {
               if( node.value.id() == newick_grammar::INNER_NODE ) {
                  // set inner node values
                  if( node.children.size() > 1 ) {
                     set_values( nodeId, ++node.children.begin(), node.children.end() );
                  }
                                                                                
                  // iterate through list children
                  tree_iterator_t i;
                  tree_container_t &childList = node.children[0].children;
                  for( i = childList.begin(); i != childList.end(); i++ ) {
                     vertex_t newNodeId = add_child( nodeId, tree.g );
                     assign_tree( newNodeId, *i );
                  }
               }
                                                                                
               if( node.value.id() == newick_grammar::LEAF_NODE ) {
                  set_values( nodeId, node.children.begin(), node.children.end() );
               }
            }


            // templated class constructor
            newick_reader_t( PTree &newTree ) : tree( newTree ) {}

            // clone this class
            base *clone() const {
               return new newick_reader_t<PTree>(tree);
            }
 
            // Read seqeunce from stream into sequence object
            void read( std::istream &input ) const {
               using namespace boost::spirit;
               using namespace boost;
               using namespace std;

               // load up until the first ';'
               string newicktree;
               getline( input, newicktree, ';' );
              if( input.fail() ) return;
               newicktree += ';';
               
               // parse the newick tree
               newick_grammar newick;
               tree_parse_info<> info = ast_parse(newicktree.c_str(), newick);

               // check if the parse was not successfull set error flag and exit.
               if( !info.match ) {
                  input.setstate( ios::failbit );
                  return;
               }
              
               // add root node implied by newick format
               tree.root = add_vertex(tree.g);

               // parse tree
               assign_tree( tree.root, info.trees[0] );
            }
         };

         // pointer to templated function
         base *bp;

      public:
         // non-templated class constructor, uses type deduction to avoid
         // gunk when calling 
         template<typename PTree>
         newick_reader( PTree &t ) {
            bp = new newick_reader_t<PTree>( t );
         }
         
         // copy constructor
         newick_reader( newick_reader const &other ) {
            bp = other.bp->clone();
         }

         // copy operator
         void operator=( newick_reader const &other ) {
            delete bp;
            bp = other.bp;
         }

         // destructor
         ~newick_reader() {
            delete bp;
         }
    
         friend std::istream &operator>>( std::istream &in, newick_reader const &fobj ) {
            fobj.bp->read( in );
            return in;
         }
   }; // end class newick_format
}

#endif
