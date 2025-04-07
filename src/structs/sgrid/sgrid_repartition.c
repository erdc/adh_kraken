#include "adh.h"
int PART_OPTION = 1;
double BALANCE_RATIO = 1.05;
double IMBALANCE_RATIO = 0.025;
int stratflag = SCOTCH_STRATDEFAULT;
//int stratflag = SCOTCH_STRATBALANCE;
//int stratflag = SCOTCH_STRATQUALITY;
//int stratflag = SCOTCH_STRATSPEED;
//int stratflag = SCOTCH_STRATSAFETY;
//int stratflag = SCOTCH_STRATREMAP;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Redistributes nodes to minimize bandwidth using PTSCOTCH library
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
int sgrid_repartition(SGRID *g) {
    int ierr = 0;
#ifdef _MESSG
    int i, icol, ierr_code = UNSET_INT;
    int npes = g->smpi->npes;
    int myid = g->smpi->myid;
    
    
    //Column flag means the grid contains columns that can not be partitioned
    // the dumb partition will respect this
    bool column_flag = FALSE;
    if (g->type == COLUMNAR) column_flag = TRUE;

    
    // Count the number of owned nodes
    int local_nnode = g->my_nnodes;
    if (column_flag) {
    	//this may need to be changed, not quite surface anymore
        local_nnode = g->my_nnodes_sur;
    }
    
    //requires node to node graph
    //may not be sufficient to just check null here, need
    //to insure it is up to date
    if (g->n2n_ind_ptr == NULL || g->n2n_edge_tab == NULL){
        sgrid_create_node_to_node_graph(g);
    }

    // flag to indicate graph weighting
    //Dont really need for SCOTCH

    // 0 for no weighting
    // 1 for edge weighting only
    // 2 for node weighting only
    // 3 for both node and edge weights
    int wgtflag = 0;
    if (column_flag) wgtflag = 2;//    // allocate the memory for the global node numbers
    int *nnode_pe = (int *) tl_alloc(sizeof(int), npes);
    int *vtxdist = (int *) tl_alloc(sizeof(int), npes + 1); // the number of nodes belonging to each pe
    
    ierr_code = MPI_Allgather(&(local_nnode), 1, MPI_INT, nnode_pe, 1, MPI_INT, g->smpi->ADH_COMM);
    if (ierr_code != MPI_SUCCESS) messg_err(ierr_code);
    
    // gather number of actual nodes on each processor to all processors
    vtxdist[0] = 0;
    for (i = 1; i <= npes; i++){
        vtxdist[i] = vtxdist[i - 1] + nnode_pe[i - 1];
    }
    
    // graph node weights if columnar, need to figure out later
    int *vwgt = NULL;
//    if (column_flag) {
//        vwgt = (int *) tl_alloc(sizeof(int), local_nnode);
//        weight_graph_nodes(g, local_nnode, vwgt);
//    }
    
    int ncon=1; //number of weights per vertex, don't think we need multi-constraint ever
    double *tpwgts = (double *) tl_alloc(sizeof(double), npes);  // fraction of vertex weight to each part ncon x npart array
    double *ubvec  = (double *) tl_alloc(sizeof(double), ncon);  // imbalance tolerance for each ncon weight
    ubvec[0] = BALANCE_RATIO; //reccomended value, could be option instead of hard coded
    for (i = 0; i < npes; i++) tpwgts[i] = 1.0/npes;//    // initialize the hash table
    

    // *xadj is exactly grid->n2n_ind_ptr
    // except we dont want the ghost node info
    // no need to copy though, it will know to stop because vtxdist array
    // maybe copy over if int gives compiler warning
    int *xadj = g->n2n_ind_ptr;

    
    // load the connections in the adjacency list array
	// adjacency already in grid->n2n_edgetab 
	// BUT need glovbal ids
	// will ommit weights for now
	// previously was physical distance but not sure that would help bandwidth
	int nedge_total = g->n2n_ind_ptr[local_nnode];
	int *adjncy = (int *) tl_alloc(sizeof(int), nedge_total); // the nodal connection table
    int *adjwgt = NULL;//(int *) tl_alloc(sizeof(int), nedge_total); // edge weights for the graph
	for (i=0; i<nedge_total; i++){
		adjncy[i] = g->node[g->n2n_edge_tab[i]].gid;
	}
    
    // partitioning
    int nparts = npes;
    int *part = (int *) tl_alloc(sizeof(int), local_nnode); // the partition returned from metis/scotch
    

    int edgecut;                       // number of edge cuts (returned)
    int numflag = 0;                   // flag to indicate numbering scheme
    int options[4] = { 0, 2, 150, 0 }; // the options for metis
    
    //parmetis call
    // maybe DM has options or something
    if (PART_OPTION == 0){
#ifdef _PARMETIS
   		assert(sizeof(real_t) == sizeof(double));
    	ierr = ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, &edgecut, part, &g->smpi->ADH_COMM);
#endif
    }else if (PART_OPTION == 1){
#ifdef _SCOTCH
    	ierr = scotch_partkway(local_nnode, xadj, adjncy, vwgt, adjwgt, stratflag, numflag, npes, IMBALANCE_RATIO, part, &g->smpi->ADH_COMM);
#endif
    }else{
    	assert(ierr==-1);
    }
    
    
    
    // cast the metis partition into part
//    int k=0, nd=UNSET_INT;
//    ID_LIST_ITEM *ptr;
//    for (i = 0; i < local_nnode; i++) {
//        if(column_flag){
//            ptr = g->vertical_list[i];
//            g->smpi->surface_partition_info[i] = (int)metis_part[i];
//            while (ptr->next != NULL) {
//                nd=ptr->id;
//                g->smpi->partition_info[nd] = (int)metis_part[i];
//                ptr = ptr->next;
//            }
//            if( g->smpi->surface_partition_info[i] != g->smpi->myid) k++;
//        }
//        else{
//            g->smpi->partition_info[i] = (int)metis_part[i];
//            if( g->smpi->partition_info[i] != g->smpi->myid) {k++;}
//        }
//    }

    for (int i = 0 ; i<local_nnode ; i++){
            printf("Rank %d part array [%d] = %d\n",myid, i, part[i]);
    }
   
    
    part = (int *) tl_free(sizeof(int), local_nnode, part);
    nnode_pe =     (int *) tl_free(sizeof(int), npes, nnode_pe);
    vtxdist =    (int *) tl_free(sizeof(int), npes + 1, vtxdist);
    adjncy =     (int *) tl_free(sizeof(int), nedge_total, adjncy);
    tpwgts =    (double *) tl_free(sizeof(double), npes, tpwgts);
    ubvec =     (double *) tl_free(sizeof(double), ncon, ubvec);
    
    if (vwgt != NULL) vwgt = (int *) tl_free(sizeof(int), local_nnode, vwgt);
#endif
    return ierr;
}
