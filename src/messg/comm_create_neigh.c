#include "adh.h"
static int DEBUG = 1;
static int MAX_NEIGH = 3; //guess for max neighboring PEs a single node could have
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
    SMPI *smpi = grid->smpi;
    //requires node to node graph
    if (grid->n2n_ind_ptr == NULL || grid->n2n_edge_tab == NULL){
        sgrid_create_node_to_node_graph(grid);
    }

    // in is easy, just resident pes of each node
    MPI_Info info = MPI_INFO_NULL;

    int reorder = 0; //no reordering of ranks
    ierr = -1;
    int indegree = 1;
    int is_unique;
    int current_node_pe;



    int myid = smpi->myid;
    int npes = smpi->npes;
    int my_nnodes = grid->my_nnodes;


    smpi->sources = tl_alloc(sizeof(int), npes);
    smpi->recv_ind = grid->my_nnodes; //starting index to recieve data
    int* source_weights_temp = tl_alloc(sizeof(int), npes);
    sarray_init_int(source_weights_temp, npes);
    //first ghost will be unique
    //KEY: ASSUMES NODES ARE ORDERED IN THAT NON-OWNED
    //NODES ARE AFTER FIRST MY_NNODES
    smpi->sources[0] = grid->node[grid->my_nnodes].resident_pe;
    //gets unique source ranks
    for (int i = grid->my_nnodes; i < grid->nnodes; i++) {
        is_unique = 1;
        current_node_pe = grid->node[i].resident_pe;
        //add to source weights
        source_weights_temp[current_node_pe]++;
        for (int j = 0; j <  indegree; j++) {
            if ( current_node_pe == grid->smpi->sources[j] ) {
                is_unique = 0; // We found a duplicate
                break;
            }
        }
        if (is_unique) {
            smpi->sources[indegree] = grid->node[i].resident_pe;
            indegree++;
        }
    }

    smpi->indegree = indegree;

    //trims vector to appropriate size 
    //should retain values? check with Corey
    smpi->sources = tl_realloc(sizeof(int), indegree, npes, grid->smpi->sources);

    //for code to work, we need sources and dests to be ordered
    assert(sarray_isordered_int(smpi->sources, indegree));

    //trims source weights to appropriate size
    smpi->source_weights = tl_alloc(sizeof(int), indegree);
    //use source weights to construct source_displs
    smpi->source_displs = tl_alloc(sizeof(int), indegree);
    smpi->source_displs[0] = 0;
    //fill in vals by indexing thorugh sources array
    for (int i = 0; i<indegree-1 ; i++){
        smpi->source_weights[i] = source_weights_temp[grid->smpi->sources[i]];
        smpi->source_displs[i+1] = smpi->source_displs[i] + smpi->source_weights[i];
    }
    smpi->source_weights[indegree-1] = source_weights_temp[grid->smpi->sources[indegree-1]];
    smpi->nrecv_neigh =  smpi->source_displs[indegree-1] + smpi->source_weights[indegree-1];
    //free memory
    source_weights_temp = tl_free(sizeof(int), npes, source_weights_temp);

    // Now need info on outgoing messages (which nodes process owns but is ghost on other process)
    // isnt it symmetric? but weights just different
    // Am i stupid or isnt it always case that if you are recieving ghost data from one PE
    // then that means you share node on an element and so you will send ur other nodes over?
    // borrowing code from comm_set_keys.c

    //loop through elems, see if any nodes belong to other PEs, this means that
    // any locally owned nodes on that element will be ghosts on other process

    //This is an element based construction similar
    // to the graph_reorder.c and the slin_sys_init_sparsity_mono
    // we could also do an nneighbor * nnode loop like in
    // comm_set_keys.c , not sure which is best


    //use n2n graph info, look at only residential nodes and who they are connected t
    int n_connections;
    int nd_start;
    int neighbor_pe;
    //use max number of edges to allocate temporary array
    int maxedges = sarray_max_distance_int(grid->n2n_ind_ptr, my_nnodes+1);
    int dest_pes[maxedges];

    int **temp_edgetab; //edge table to store which PEs residential node may send to
    //allocate based on assumption of MAX_NEIGH, will realloc if not true
    temp_edgetab = (int**) tl_alloc(sizeof(int*), my_nnodes);
    for(int j=0;j<my_nnodes;j++){
        temp_edgetab[j] = (int*) tl_alloc(sizeof(int), MAX_NEIGH);
        for(int k=0;k<MAX_NEIGH;k++){
            //initializing with INT_MAX to make sorting easier
            temp_edgetab[j][k]=INT_MAX;
        }
    }


    int* dest_weights_temp = tl_alloc(sizeof(int), npes);
    sarray_init_int(dest_weights_temp,npes);


    //also create array of sizes for each, may be necessary to store
    int *temp_nedges = (int*) tl_alloc(sizeof(int), my_nnodes);
    sarray_init_int(temp_nedges, my_nnodes); // start at 0

    int prev_pe;
    int ctr;
    int n_send;
    
    for (int i = 0; i < grid->my_nnodes; i++){
        n_connections = grid->n2n_ind_ptr[i+1] - grid->n2n_ind_ptr[i];
        nd_start = grid->n2n_ind_ptr[i];
        prev_pe = UNSET_INT;
        ctr=0;
        sarray_init_value_int(dest_pes, maxedges, INT_MAX);
        //loop over edges and find which unique PEs it will be sending to
        for(int j = 0; j<n_connections ; j++){
            //current connection
             neighbor_pe =  grid->node[grid->n2n_edge_tab[nd_start+j]].resident_pe;
             //if current connection is not residential then will need to send info to them
             if( neighbor_pe != myid && neighbor_pe != prev_pe){

                dest_pes[ctr] = neighbor_pe;
                ctr++;
                prev_pe = neighbor_pe;
            }
        }

        //we have other PEs dest_pes stored with possible duplicates, need to take unique
        //if all neighbors are residential skip
        if( ctr!= 0){
            //sort first or else unique wont work
            qsort(dest_pes, ctr, sizeof(int), compare_ints);
            n_send = sarray_unique_int(dest_pes, ctr);
            temp_nedges[i] = n_send;

            //if n_send is > MAX_NEIGH then realloc
            //should happen very seldom if ever
            if (n_send > MAX_NEIGH){
                temp_edgetab[i] = (int*) tl_realloc(sizeof(int), n_send, MAX_NEIGH, temp_edgetab[i]);
            }
            //store the dest_pes into the temp_edgetab
            for (int j = 0; j < n_send; j++){
                temp_edgetab[i][j] = dest_pes[j];
                //add to dest_weights_temp
                dest_weights_temp[dest_pes[j]] += 1;
            }

        }


    }


    //analyze dest_weights_temp to get outdegree, dest_weights, dest
    smpi->outdegree =  sarray_num_nonzero_int(dest_weights_temp, npes);

    //for ghost communication, neighbors should be symmetric
    assert(smpi->outdegree == smpi->indegree );

    //allocate and fill out smpi->dest
    smpi->dest = tl_alloc(sizeof(int), smpi->outdegree);
    assert(sarray_get_indeces_nonzero_int(smpi->dest, dest_weights_temp, npes) == smpi->outdegree);

    //trims source weights to appropriate size
    smpi->dest_weights = tl_alloc(sizeof(int), smpi->outdegree);
    //use source weights to construct source_displs
    smpi->dest_displs = tl_alloc(sizeof(int), smpi->outdegree);
    smpi->dest_displs[0] = 0;
    //fill in vals by indexing thorugh sources array
    for (int i = 0; i < (smpi->outdegree-1) ; i++){
        smpi->dest_weights[i] = dest_weights_temp[smpi->dest[i]];
        smpi->dest_displs[i+1] = smpi->dest_displs[i] + smpi->dest_weights[i];
    }
    smpi->dest_weights[smpi->outdegree-1] = dest_weights_temp[smpi->dest[smpi->outdegree-1]];
    


    //free memory
    dest_weights_temp = tl_free(sizeof(int), npes, dest_weights_temp);

    //use temp_edgetab to fill out the actual table we want
    int n_indices = smpi->dest_displs[smpi->outdegree-1] + smpi->dest_weights[smpi->outdegree-1];
    smpi->nsend_neigh = n_indices;
    smpi->dest_indices = (int*) tl_alloc(sizeof(int), n_indices);

    int *temp_map = (int *) tl_alloc(sizeof(int), npes);
    sarray_init_value_int(temp_map, npes, UNSET_INT);
    for(int i = 0; i < smpi->outdegree ; i++){
        temp_map[ smpi->dest[i] ] = i;
    }

    int neighbor_no;
    int *temp_count = tl_alloc(sizeof(int), smpi->outdegree);
    sarray_init_int(temp_count, smpi->outdegree);


    for(int i = 0; i < my_nnodes; i++){
        for(int j = 0; j < temp_nedges[i]; j++){
            neighbor_no = temp_map[temp_edgetab[i][j]];
            smpi->dest_indices[ smpi->dest_displs[neighbor_no] + temp_count[neighbor_no] ] = i;
            temp_count[neighbor_no]++;
        }
    }


    //for(int i = 0; i<n_indices; i++ ){
    //    printf("Rank [%d] dest_indices [%d] = %d \n", myid, i, smpi->dest_indices[i]);
    //}
    // do we need temp_edgetab? for now i guess not. We have structure in snode to stor
    // for now don't see need


    //clear out temporary variables
    for(int j=0;j<my_nnodes;j++){
        temp_edgetab[j] = (int*) tl_free(sizeof(int), MAX(temp_nedges[j],MAX_NEIGH), temp_edgetab[j]);
    }
    temp_edgetab = (int**) tl_free(sizeof(int*), my_nnodes,temp_edgetab);
    temp_nedges = (int*) tl_free(sizeof(int), my_nnodes, temp_nedges);
    temp_map = (int *) tl_free(sizeof(int), npes, temp_map);
    temp_count = (int*) tl_free(sizeof(int), smpi->outdegree, temp_count);





    //could optionally include source weights or assume all equal with MPI_UNWEIGHTED
    ierr =  MPI_Dist_graph_create_adjacent(smpi->ADH_COMM,
                                           smpi->indegree, smpi->sources,
                                           smpi->source_weights,
                                           smpi->outdegree, smpi->dest,
                                           smpi->dest_weights,
                                           info, reorder, &(smpi->ADH_NEIGH));




    if (DEBUG){
        assert(grid->smpi->ADH_NEIGH);

        int neighbors;
        int *in_neighbors = NULL, *out_neighbors = NULL;
        int *in_weights = NULL, *out_weights = NULL;

        MPI_Dist_graph_neighbors_count(smpi->ADH_NEIGH, &(smpi->indegree), &(smpi->outdegree), &neighbors);

        in_neighbors = (int*) tl_alloc(sizeof(int), smpi->indegree);
        in_weights = (int*) tl_alloc(sizeof(int), smpi->indegree);
        out_neighbors = (int*) tl_alloc(sizeof(int), smpi->outdegree);
        out_weights = (int*) tl_alloc(sizeof(int), smpi->outdegree);

        MPI_Dist_graph_neighbors(smpi->ADH_NEIGH, smpi->indegree, in_neighbors, in_weights,
                                 smpi->outdegree, out_neighbors, out_weights);

        printf("Rank %d: In-degree = %d, Out-degree = %d\n", myid, smpi->indegree, smpi->outdegree);
        printf("Rank %d: In-neighbors: ", myid);
        for (int i = 0; i < smpi->indegree; i++) {
            printf("(%d, weight=%d) ", in_neighbors[i], in_weights[i]);
        }
        printf("\n");
        printf("Rank %d: Out-neighbors: ", myid);
        for (int i = 0; i < smpi->outdegree; i++) {
            printf("(%d, weight=%d) ", out_neighbors[i], out_weights[i]);
        }
        printf("\n");
    }

#endif
    return ierr;

}
