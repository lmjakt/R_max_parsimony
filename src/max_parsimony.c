#include <R.h>
#include <Rinternals.h>
#include "tree.h"

// tree_r        : a matrix with two columns: parent, child
//                 that describes connections within the tree
// tree_prop_r   : the number of nodes and leaf_nodes
// sub_matrix_r  : a substitution matrix.
// alphabet_r    : al_offset and al_size used for the substitution matrix
// leaf_states_r : a character vector, where each entry defines the
//                 state of a leaf_node.
//  Note that how missing values are handled depends on the substitution
//  tree and 
SEXP sankoff(SEXP tree_r, SEXP tree_props_r, SEXP sub_matrix_r, SEXP alphabet_r,
	     SEXP leaf_states_r)
{
  if(TYPEOF(tree_r) != INTSXP || TYPEOF(tree_props_r) != INTSXP ||
     TYPEOF(sub_matrix_r) != INTSXP || TYPEOF(alphabet_r) != INTSXP
     || TYPEOF(leaf_states_r) != STRSXP){
    error("arguments should all be integers except for leaf_states which should be character");

    SEXP tree_dims_r = getAttrib( tree_r, R_DimSymbol );
    int *tree_dims = INTEGER(tree_dims_r);
    if(tree_dims[2] != 2)
      error("tree should have two columns");

    if(length(tree_props_r) != 2)
      error("tree_props_r should have a length of 2");
    int *tree_props = INTEGER(tree_props_r);
    int nodes_n = tree_props[0];
    int leaf_n = tree_props[1];
    if(leaf_n >= nodes_n || leaf_n < 2)
      error("there must be at least two leaves and as many non-leaf nodes");

    SEXP sub_matrix_dims_r = getAttrib( sub_matrix_r, R_DimSymbol );
    if(length(sub_matrix_dims_r) != 2)
      error("the substitution matrix must be a two dimensional matrix");
    int *sub_matrix_dims = INTEGER(sub_matrix_dims_r);

    if(length(alphabet_r) != 2)
      error("The alphabet description must contain two words");
    int al_offset = INTEGER(alphabet_r)[0];
    int al_size = INTEGER(alphabet_r)[1];
   

}
