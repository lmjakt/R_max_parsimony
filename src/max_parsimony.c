#include <R.h>
#include <Rinternals.h>
#include "tree.h"

// tree_r        : a matrix with two columns: parent, child
//                 that describes connections within the tree
// tree_prop_r   : the number of nodes and leaf_nodes
// sub_matrix_r  : a substitution matrix.
// al_dims       : al_offset and al_size used for the substitution matrix
// leaf_states_r : a character vector, where each entry defines the
//                 state of a leaf_node.
//  Note that how missing values are handled depends on the substitution
//  tree and 
SEXP sankoff(SEXP tree_r, SEXP tree_props_r, SEXP sub_matrix_r, SEXP al_dims,
	     SEXP leaf_states_r)
{


}
