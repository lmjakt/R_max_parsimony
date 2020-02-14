#include <stdlib.h>
#include <string.h>
#include <Rinternals.h>  // for Rprintf()
#include "tree.h"

// this is the biggest value that can be multiplied by 2 and still give
// a positive number. We use signed int as R does not have unsigned
// integers.
const int max_length = (int)(((unsigned int)~0) >> 2);

// returns the number of OK characters. These should always enumerate to 
// dim_n.
int check_state(const char *desc, int al_offset, int al_size){
  int i = 0;
  while(desc[i] != 0){
    unsigned char c = (unsigned char)desc[i];
    if( c < al_offset || c >= (al_offset + al_size) )
      return(i);
    ++i;
  }
  return(i);
}

// confirm that the information provided is correct... 
void check_tree(struct h_tree *tree){
  tree->is_good = true;
  for(int i=0; i < tree->leaf_n; ++i){
    int c = check_state( tree->leaf_states[i], tree->al_offset, tree->al_size );
    if(c != tree->dim_n){
      tree->is_good = false;
      return;
    }
  }
  // also check that there are leaf_n leaf nodes and that
  // these are numbered from 1:leaf_n in the edge_child array..
  // do this by counting the number of children
  unsigned int *child_count = malloc(sizeof(unsigned int) * tree->node_n);
  unsigned int *parent_count = malloc(sizeof(unsigned int) * tree->node_n);
  memset( child_count, 0, sizeof(unsigned int) * tree->node_n );
  memset( parent_count, 0, sizeof(unsigned int) * tree->node_n );
  // here parent_count means the number of times an index appears in the
  // the parent array. That is equivalent to the number of children a
  // given parent has. 
  for(int i=0; i < tree->edge_n; i++){
    int p_i = tree->edge_parent[i] - 1;
    int c_i = tree->edge_child[i] - 1;
    if(p_i >= 0 && c_i >= 0 && p_i < tree->node_n && c_i < tree->node_n){
      parent_count[p_i]++;
      child_count[c_i]++;
    }else{
      tree->is_good = false;
      break;  // don't return as we have to free memory
    }
  }
  // then go through and count the number of nodes which do not have
  // children and make sure that these have appropriate indices
  // We can also count the number of parents
  int n = 0;
  for(int i=0; i < tree->node_n; ++i){
    if( parent_count[i] == 0 ){
      n++;
      if( i >= tree->leaf_n )
	tree->is_good = false;
    }
  }
  // and finally make sure that the n is equal to the specified leaf node number
  if( n != tree->leaf_n )
    tree->is_good = false;
  free(child_count);
  free(parent_count);
}

struct h_tree make_tree(int *edge_child, int *edge_parent, int edge_n, int node_n, int leaf_n,
			int dim_n, const char **leaf_states, int al_offset, int al_size){
  struct h_tree tree;
  tree.node_n = node_n;
  tree.leaf_n = leaf_n;
  tree.edge_n = edge_n;
  tree.edge_child = edge_child;
  tree.edge_parent = edge_parent;
  // the leaf_states are in order.
  tree.dim_n = dim_n;
  tree.leaf_states = leaf_states;
  tree.al_offset = al_offset;
  tree.al_size = al_size;
  //
  check_tree( &tree );
  return(tree);
}

// Return an array of nodes..
// Where the leaf nodes have set the
struct ht_node* make_nodes(struct h_tree *tree, int *sub_matrix, int *root_i){
  struct ht_node *nodes = malloc( sizeof(struct ht_node) * tree->node_n );
  memset( nodes, 0, sizeof(struct ht_node) * tree->node_n );
  // then we go through the parent and child arrays.. and assign.
  for(int i=0; i < tree->edge_n; ++i){
    int p_i = tree->edge_parent[i] - 1;
    int c_i = tree->edge_child[i] - 1;
    if(p_i >= 0 && c_i >= 0 && p_i < tree->node_n && p_i < tree->node_n){
      if( nodes[p_i].edge_n < 3 ){
	nodes[p_i].edges[ nodes[p_i].edge_n ] = &nodes[c_i];
	nodes[p_i].edges_i[ nodes[p_i].edge_n ] = c_i;
	nodes[p_i].is_child[ nodes[p_i].edge_n ] = true;
	nodes[p_i].edge_n++;
      }
      if( nodes[c_i].edge_n < 3){
	nodes[c_i].edges[ nodes[c_i].edge_n ] = &nodes[p_i];
	nodes[c_i].edges_i[ nodes[c_i].edge_n ] = p_i;
	nodes[c_i].is_child[ nodes[c_i].edge_n ] = false;
	nodes[c_i].edge_n++;
      }
    }
  }
  // 
  // Set up the leaf states;
  *root_i = -1;
  for(int i=0; i < tree->node_n; ++i){
    // check if we have a root.. (root has no parent)
    int p_count = 0;
    for(int j=0; j < nodes[i].edge_n; ++j)
      p_count += (nodes[i].is_child[j] == true ? 0 : 1);
    if(p_count == 0)
      *root_i = i;
    nodes[i].tree_lengths = malloc( sizeof(unsigned int) * tree->al_size * tree->dim_n );
    nodes[i].length_determined = false;
    if(i < tree->leaf_n){  // we know the state..
      // Default to setting all the values to a rather large number,
      // that can just about be multiplied by 2... 
      for(int j=0; j < (tree->al_size * tree->dim_n); ++j)
	nodes[i].tree_lengths[j] = max_length;
      for(int j=0; j < tree->dim_n; ++j){
	int o = tree->leaf_states[i][j] - tree->al_offset;
	nodes[i].tree_lengths[ j * tree->al_size + o  ] = 0;
      }
      nodes[i].length_determined = true;
    }
  }
  // return as an array, since the root has been set we can leave it to the caller to handle everything here..
  return( nodes );
}

void ht_nodes_free(struct ht_node *nodes, int l){
  for(int i=0; i < l; ++i){
    free(nodes[i].tree_lengths);
  }
  free(nodes);
}


void sankoff_set_lengths( struct ht_node *node, int *sub_matrix, int al_offset, int al_size, int dim_n ){
  if( node->length_determined )
    return;
  
  // at some point we want to merge the information from two children..
  struct ht_node *children[3];
  memset( children, 0, sizeof( struct ht_node* ) * 3 );
  int child_count = 0;
  for( int i=0; i < 3; ++i ){
    if( node->edges[i] && node->is_child[i] ){
      sankoff_set_lengths( node->edges[i], sub_matrix, al_offset, al_size, dim_n );
      if(node->edges[i]->length_determined){
	children[ child_count ] = node->edges[i];
	child_count++;
      }
    }
  }
  // if we have no children we can not have any state;
  if(child_count == 0){
    node->length_determined = false;
    return;
  }
  // if we have only one child this is not a reasonable tree, but nevertheless we should simply
  // have exactly the same state as the child.
  if(child_count == 1){
    memcpy( node->tree_lengths, children[0]->tree_lengths, sizeof(int) * al_size * dim_n );
    node->length_determined = true;
  }

  // In order to handle an arbitrary number of children we define
  // an array of pointers to the children length states and
  // an array ints where we store the minimal length values.
  int **c_l = malloc( sizeof(int*) * child_count);
  int *ll = malloc(sizeof(int) * child_count);

  // And merge the tree lengths
  for(int i=0; i < dim_n; ++i){
    for(int j=0; j < child_count; ++j)
      c_l[j] = children[j]->tree_lengths + i * al_size;
    int *l3 = node->tree_lengths + i * al_size;
    // determine l3 from the values in c_l
    for(int j=0; j < al_size; ++j){  // j -> position in node
      for(int k=0; k < child_count; ++k)
	ll[k] = max_length;
      for(int k=0; k < al_size; ++k){ // k -> position in children
	int cost = sub_matrix[ j * al_size + k ];
	for(int m=0; m < child_count; ++m){
	  if( ll[m] > c_l[m][k] + cost )
	    ll[m] = c_l[m][k] + cost;
	}
      }
      l3[j] = 0;
      for(int k=0; k < child_count; ++k)
	l3[j] += ll[k];
    }
  }
  node->length_determined = true;
  free( c_l );
  free( ll );
  // We can then go up this tree in the other direction. And make a prediction.
  return;
}
