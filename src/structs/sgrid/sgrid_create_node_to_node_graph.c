#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes node-to-node connectivity
 * 			   Useful in CG sparsity pattern, graph reordering with SCOTCH, or
 *             establishing neighborhood communicators
 * 			   will fill out grid->n2n_* stuff
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in,out] pgrid (SGRID *)  pointer to an AdH grid
 * @returns int - 0 if succesful
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int sgrid_create_node_to_node_graph(SGRID *grid){

	int i,j,k;
	int err = 0;
	int n_connections;
	int nelems3d = grid->nelems3d;
	int nelems2d = grid->nelems2d;
	int nelems1d = grid->nelems1d;
	int current_node, other_node;
	int count=0;
	int nnodes = grid->nnodes;
	int nnodes_on_elem;
	int *nodes; //for aliasing
	int Nedges=0;
	int *n_con;
	int *n_con_no_dup;
	int **temp_edgetab;
	//quickly compute number of edges in the mesh
	//the nodes will not be as trivial, need to know how many elements each node is connected to
	//loop through elements and add up
	//need nodal array to store this

	//need two nodal arrays
	//these will store the number of nodes each node is connected to
	n_con = (int *) tl_alloc(sizeof(int), nnodes);
    n_con_no_dup = (int *) tl_alloc(sizeof(int), nnodes);
    grid->n2n_ind_ptr = (int *) tl_alloc(sizeof(int),nnodes+1);
    sarray_init_int(n_con, nnodes);
    sarray_init_int(n_con_no_dup, nnodes);
    //First set of loops is solely to establish how many nodes are connected to each node
    //could avoid first pass by guessing max connections per node
    for (i=0;i<nelems3d;i++){
    	//nnodes on the element
        nnodes_on_elem = grid->elem3d[i].nnodes;
        //Could get stuff for node weights from physics mat
       	//for this element, find nodes
        nodes = grid->elem3d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
            	//other_node = nodes[k];
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
            	//other_node = nodes[k];
            	n_con[current_node]+=1;
            }
        }
    }
    for (i=0;i<nelems2d;i++){
    	//nnodes on the element
        nnodes_on_elem = grid->elem2d[i].nnodes;
        //Could get stuff for node weights from physics mat
       	//for this element, find nodes
        nodes = grid->elem2d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
            	//other_node = nodes[k];
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
            	//other_node = nodes[k];
            	n_con[current_node]+=1;
            }
        }
    }
    for (i=0;i<nelems1d;i++){
    	//nnodes on the element
        nnodes_on_elem = grid->elem1d[i].nnodes;
        //Could get stuff for node weights from physics mat
       	//for this element, find nodes
        nodes = grid->elem1d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
            	//other_node = nodes[k];
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
            	//other_node = nodes[k];
            	n_con[current_node]+=1;
            }
        }
    }
    //use nnz_rows to dynamically allocate
    //int temp_cols_diag[nrows][nCon3d];
    temp_edgetab = (int**) tl_alloc(sizeof(int*), nnodes);
    for(j=0;j<nnodes;j++){
        temp_edgetab[j] = (int*) tl_alloc(sizeof(int), n_con[j]);
        for(k=0;k<n_con[j];k++){
            temp_edgetab[j][k]=INT_MAX;
        }
    }
    //Seems redundant but must reuse as indexing
    sarray_init_int(n_con, nnodes);
    //Now fill in edgetab, will contain duplicates
    for (i=0;i<nelems3d;i++){
    	//nnodes on the element
        nnodes_on_elem = grid->elem3d[i].nnodes;
        //Could get stuff for node weights from physics mat
       	//for this element, find nodes
        nodes = grid->elem3d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
            	other_node = nodes[k];
                temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
            	other_node = nodes[k];
            	temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
        }
    }
    for (i=0;i<nelems2d;i++){
    	//nnodes on the element
        nnodes_on_elem = grid->elem2d[i].nnodes;
        //Could get stuff for node weights from physics mat
       	//for this element, find nodes
        nodes = grid->elem2d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
            	other_node = nodes[k];
                temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
            	other_node = nodes[k];
            	temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
        }
    }
    for (i=0;i<nelems1d;i++){
    	//nnodes on the element
        nnodes_on_elem = grid->elem1d[i].nnodes;
        //Could get stuff for node weights from physics mat
       	//for this element, find nodes
        nodes = grid->elem1d[i].nodes;
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (j=0;j<nnodes_on_elem;j++){
            //each node inside an element, add this to the count
            current_node = nodes[j];
            //loop through each node that is not the node itself
            //nodes before current node
            for (k=0;k<j;k++){
            	other_node = nodes[k];
                temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
            for(k=j+1;k<nnodes_on_elem;k++){
            	other_node = nodes[k];
            	temp_edgetab[current_node][n_con[current_node]]=other_node;
                n_con[current_node]+=1;
            }
        }
    }
    //sort this and eliminate duplicates
    for (i=0;i<nnodes;i++){
        // sort the column indices (j-entries)
        //use stdlib.h qsort
        qsort(temp_edgetab[i], n_con[i], sizeof(int), compare_ints);
        //this should hopefully remove duplicates?
        n_connections = sarray_unique_int(temp_edgetab[i], n_con[i]);
        //overwrite nnz row with sarray_unique_int?
        n_con_no_dup[i] = n_connections;
        Nedges+=n_connections;
    }
    //allocate adjacency data
    grid->n2n_edge_tab = (int *) tl_alloc(sizeof(int), Nedges);
    int *nodal_edges; //alias

    //now use info to fill in the edgetan;e
    for(i=0;i<nnodes;i++){
        //printf("filling in index ptr and column entries, row %d\n",i);
        grid->n2n_ind_ptr[i] = count;
        //printf("filling in index ptr and column entries, row %d\n",i);
        nodal_edges = temp_edgetab[i];
        //think about how to do this since each temp_rows may be different size
        for(j=0;j<n_con_no_dup[i];j++){
            grid->n2n_edge_tab[count] = nodal_edges[j];
            count++;
        }
    }
    //also last entry
    grid->n2n_ind_ptr[nnodes] = count;
    assert(count == Nedges);

    //freeing memory
    n_con_no_dup = (int *) tl_free(sizeof(int), nnodes, n_con_no_dup);
    for(j=0;j<nnodes;j++){
        temp_edgetab[j] = (int*) tl_free(sizeof(int), n_con[j],temp_edgetab[j]);
    }
    temp_edgetab = (int**) tl_free(sizeof(int*), nnodes,temp_edgetab);
    n_con = (int *) tl_free(sizeof(int), nnodes, n_con);

    return err;
}
