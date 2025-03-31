#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Establish neighborhood communicators
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
int comm_create_neighborhood(SGRID *grid){
    int ierr = 0;
#ifdef _MESSG
    // Define the local adjacencies for each process
    // using ghost nodes

    // To establish neighborhood communicator
    // need # processes it will recieve from
    // and # processes it will send to

    // in is easy, just resident pes of each node
    MPI_Info info = MPI_INFO_NULL;
    int reorder = 0; //no reordering of ranks
    ierr = -1;
    int indegree = 1;
    int is_unique;
    int current_node_pe;
    int myid = grid->smpi->myid;
    int npes = grid->smpi->npes;
    int* sources = tl_alloc(sizeof(int), npes);
    int* source_weights_temp = tl_alloc(sizeof(int), npes);
    sarray_init_int(source_weights_temp, npes);
    //first ghost will be unique
    //KEY: ASSUMES NODES ARE ORDERED IN THAT NON-OWNED
    //NODES ARE AFTER FIRST MY_NNODES
    sources[0] = grid->node[grid->my_nnodes].resident_pe;
    //gets unique source ranks
    for (int i = grid->my_nnodes; i < grid->nnodes; i++) {
        is_unique = 1;
        current_node_pe = grid->node[i].resident_pe;
        //add to source weights
        source_weights_temp[grid->node[i].resident_pe]++;
        for (int j = 0; j <  indegree; j++) {
            if ( current_node_pe == sources[j] ) {
                is_unique = 0; // We found a duplicate
                break;
            }
        }
        if (is_unique) {
            sources[indegree] = grid->node[i].resident_pe;
            indegree++;
        }
    }

    //trims vector to appropriate size 
    //should retain values? check with Corey
    sources = tl_realloc(sizeof(int), npes, indegree, sources);

    //trims source weights to appropriate size
    int *source_weights = tl_alloc(sizeof(int), indegree);
    //fill in vals by indexing thorugh sources array
    for (int i = 0; i<indegree ; i++){
        source_weights[i] = source_weights_temp[sources[i]];
    }

    //free memory
    source_weights_temp = tl_free(sizeof(int), npes, source_weights_temp);


    //could optionally include source weights or assume all equal with MPI_UNWEIGHTED
    //for now just don't weight


    // Now need info on outgoing messages (which nodes process owns but is ghost on other process)
    // isnt it symmetric? but weights just different
    // Am i stupid or isnt it always case that if you are recieving ghost data from one PE
    // then that means you share node on an element and so you will send ur other nodes over?
    // borrowing code from comm_set_keys.c

    //loop through elems, see if any nodes belong to other PEs, this means that
    // any locally owned nodes on that element will be ghosts on other process
    int *dest_weights_temp = tl_alloc(sizeof(int), npes);
    sarray_init_int(dest_weights_temp, npes);
    

    int ie; int inode;
    
    int temp_arr[MAX_NNODE];
    int temp_ind_ptr[MAX_NNODE];
    int nnode;
    int elem_has_n_ghost_nodes;
    int elem_has_n_ghost_pes;
    int n_neighbors;
    int n_residential;



    for (ie = 0; ie < grid->nelems3d; ie++){
        nnode = grid->elem3d[ie].nnodes;
        elem_has_n_ghost_nodes = 0;
        n_neighbors = 0;
        n_residential = 0;
        //store residential PEs into temporary array
        for (inode = 0; inode <  nnode; inode++){     
            current_node_pe = grid->node[grid->elem3d[ie].nodes[inode]].resident_pe;

            //store residential nodes only
            if (current_node_pe == myid){ 
                // if inode'th node isnt residential
                // then this element contains ghost
                // and so all residential nodes will be sent to other PEs
                // coule be multiple ghosts from different PEs
                temp_ind_ptr[n_residential] = inode;
                n_residential += 1;
            }else{
                //if there are ghost nodes present, keep track of the PEs
                temp_arr[elem_has_n_ghost_nodes] = current_node_pe;
                elem_has_n_ghost_nodes += 1;
            }
        }

        //if all nodes are residential then skip
        if (elem_has_n_ghost_nodes > 0){
            // eliminate duplicates in temp_arr
            elem_has_n_ghost_pes = sarray_unique_int(temp_arr, elem_has_n_ghost_nodes);
            // use temp_arr to generate dest and dest weights
            for(inode = 0; inode < elem_has_n_ghost_pes; inode++){
                //for each unique non-residential PE
                //add # of residential nodes to dest weights
                dest_weights_temp[temp_arr[inode]] += n_residential;
            }
        }

    }

    for (ie = 0; ie < grid->nelems2d; ie++){
        nnode = grid->elem2d[ie].nnodes;
        elem_has_n_ghost_nodes = 0;
        n_neighbors = 0;
        n_residential = 0;
        //store residential PEs into temporary array
        for (inode = 0; inode <  nnode; inode++){     
            current_node_pe = grid->node[grid->elem2d[ie].nodes[inode]].resident_pe;

            //store residential nodes only
            if (current_node_pe == myid){ 
                // if inode'th node isnt residential
                // then this element contains ghost
                // and so all residential nodes will be sent to other PEs
                // coule be multiple ghosts from different PEs
                temp_ind_ptr[n_residential] = inode;
                n_residential += 1;
            }else{
                //if there are ghost nodes present, keep track of the PEs
                temp_arr[elem_has_n_ghost_nodes] = current_node_pe;
                elem_has_n_ghost_nodes += 1;
            }
        }

        //if all nodes are residential then skip
        if (elem_has_n_ghost_nodes > 0){
            // eliminate duplicates in temp_arr
            elem_has_n_ghost_pes = sarray_unique_int(temp_arr, elem_has_n_ghost_nodes);
            // use temp_arr to generate dest and dest weights
            for(inode = 0; inode < elem_has_n_ghost_pes; inode++){
                //for each unique non-residential PE
                //add # of residential nodes to dest weights
                dest_weights_temp[temp_arr[inode]] += n_residential;
            }
        }

    }

    for (ie = 0; ie < grid->nelems1d; ie++){
        nnode = grid->elem1d[ie].nnodes;
        elem_has_n_ghost_nodes = 0;
        n_neighbors = 0;
        n_residential = 0;
        //store residential PEs into temporary array
        for (inode = 0; inode <  nnode; inode++){     
            current_node_pe = grid->node[grid->elem1d[ie].nodes[inode]].resident_pe;

            //store residential nodes only
            if (current_node_pe == myid){ 
                // if inode'th node isnt residential
                // then this element contains ghost
                // and so all residential nodes will be sent to other PEs
                // coule be multiple ghosts from different PEs
                temp_ind_ptr[n_residential] = inode;
                n_residential += 1;
            }else{
                //if there are ghost nodes present, keep track of the PEs
                temp_arr[elem_has_n_ghost_nodes] = current_node_pe;
                elem_has_n_ghost_nodes += 1;
            }
        }

        //if all nodes are residential then skip
        if (elem_has_n_ghost_nodes > 0){
            // eliminate duplicates in temp_arr
            elem_has_n_ghost_pes = sarray_unique_int(temp_arr, elem_has_n_ghost_nodes);
            // use temp_arr to generate dest and dest weights
            for(inode = 0; inode < elem_has_n_ghost_pes; inode++){
                //for each unique non-residential PE
                //add # of residential nodes to dest weights
                dest_weights_temp[temp_arr[inode]] += n_residential;
            }
        }

    }

    //after all elements are traversed
    //use dest_weights_temp to fill out outdegree, dest, and destweights

    //outdegree will be number of nonzeros in dest_weights_temp
    int outdegree =  sarray_num_nonzero_int(dest_weights_temp, npes);
    //should always be same as indegree
    assert(outdegree == indegree);
    int *dest = tl_alloc(sizeof(int), outdegree);
    int *dest_weights = tl_alloc(sizeof(int), outdegree);

    sarray_get_indeces_nonzero_int(dest, dest_weights_temp, npes);
    for(int i = 0; i < outdegree; i++){
        dest_weights[i] = dest_weights_temp[dest[i]];
    }

    //free memory
    dest_weights_temp = tl_free(sizeof(int), npes, dest_weights_temp);



    ierr =  MPI_Dist_graph_create_adjacent(grid->smpi->ADH_COMM,
                                   indegree, sources,
                                   source_weights,
                                   outdegree, dest,
                                   dest_weights,
                                   info, reorder, &(grid->smpi->ADH_NEIGH));
#endif
    return ierr;

}