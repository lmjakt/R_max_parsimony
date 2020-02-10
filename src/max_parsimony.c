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
  // apart from leaf_states_r, everything should be INTSXP. leaf_states_r should be
  // STRSXP
  if(TYPEOF(tree_r) != INTSXP || TYPEOF(tree_props_r) != INTSXP || TYPEOF(sub_matrix_r) != INTSXP
     || TYPEOF(al_dims) != INTSXP || TYPEOF(leaf_states_r) != STRSXP )
    error("incorrect data type, should be: int, int, int, int, character");
  SEXP tree_dims = getAttrib(tree_r, R_DimSymbol);
  if(length(tree_dims) != 2)
    error("the tree should be a two-dimensional matrix");
  int tree_nrow = INTEGER(tree_dims)[0];
  int tree_ncol = INTEGER(tree_dims)[1];
  if(tree_ncol != 2 || tree_nrow < 3)
    error("The tree matrix should have at least 3 rows and exactly two columns");
  int *tree_edges = INTEGER(tree_r);
  int *tree_parents = tree;
  int *tree_children = tree + tree_nrow;

  // tree_props_r should have length 2
  if(length(tree_props_r) != 2)
    error("tree_props should have a length of 2");
  int node_n = INTEGER(tree_props_r)[0];
  int leaf_n = INTEGER(tree_props_r)[1];
  if(leaf_n < 2 || node_n <= leaf_n)
    error("there should be more nodes than leaves");

  // sub_matrix should be a square matrix
  SEXP sub_matrix_dims_r = getAttrib(sub_matrix_r, R_DimSymbol);
  if(length(sub_matrix_dims_r) != 2)
    error("sub_matrix should be a matrix");
  int *sub_matrix_dims = INTEGER(sub_matrix_dims_r);
  if(sub_matrix_dims[0] != sub_matrix_dims[1])
    error("sub_matrix should be square");
  int *sub_matrix = INTEGER(sub_matrix_r);

  // al dims should be an integer vector of length 2
  if(length(al_dims) != 2)
    error("al_dims should have length 2");
  int al_offset = INTEGER(al_dims)[0];
  int al_size = INTEGER(al_dims)[1];

  if(al_size != sub_matrix_dims[0])
    error("the dimensions of the submatrix should be the same as the alphabet size");

  // leaf_states_r should have the same length as the number of leaf nodes specified
  if(length(leaf_states_r) != leaf_n)
    error("the length of leaf_states_r should be the same as leaf_n");

  int dim_n = length( STRING_ELT(leaf_states_r, 0));
  const char **leaf_states = malloc(sizeof(char*) * length(leaf_states_r));
  for(int i=0; i < length(leaf_states_r); ++i){
    SEXP s = STRING_ELT(leaf_states_r, i);
    if(length(s) != dim_n){
      free(leaf_states);
      error("all leaf_states must be the same length");
    }
    leaf_states[i] = CHAR( s );
  }
  // and we should now be able to make the tree and then infer it's ancestry...

  struct h_tree tree = make_tree( tree_children, tree_parents, tree_nrow, node_n, leaf_n,
				  dim_n, leaf_states, al_offset, al_size );
  int root_i = -1;
  struct ht_node* nodes = make_nodes( &tree, sub_matrix, &root_i);
  if(root_i >= 0)
    int n = sankoff_set_lengths( nodes + root_i, sub_matrix, al_offset, al_size, dim_n );

  // That should give us a tree with distances in all the roots. We haven't yet solved
  // the issue of what to do with the root that has three children but no parents.
  // But in any case it will be necessary to somehow validat the tree before we do anything
  // else.

  // But lets see if we can compile it first.
  ht_node
  
}
