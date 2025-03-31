//#include "adh.h"
///*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
///*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
///*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
///*!
// *  \brief     Redistributes nodes to minimize bandwidth using PTSCOTCH library
// *  \author    Corey Trahan, Ph.D.
// *  \author    Mark Loveand, Ph.D.
// *  \bug       none
// *  \warning   none
// *  \copyright AdH
// *
// * @param[in,out] pgrid (SGRID *)  pointer to an AdH grid
// * @returns int - 0 if succesful
// *
// */
///*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//int sgrid_repartition_ptscotch(SGRID *g) {
//    int i, icol, ierr_code = UNSET_INT;
//    EDGE_LIST_ITEM *edge_hashtab[HASHSIZE];   // hash table for the edges
//    int npes = g->smpi->npes;
//    int myid = g->smpi->myid;
//    
//    //Mark, may need to rethink here
//    bool column_flag = FALSE;
//    if (g->type == COLUMNAR) column_flag = TRUE;
//    //assert(g->type == COLUMNAR);
//    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
//    
//    // Count the number of owned nodes
//    int local_nnode = g->my_nnodes;
//    if (column_flag) {
//        local_nnode = g->my_nnodes_sur;
//    }
//    
//    // flag to indicate graph weighting
//    idx_t wgtflag = 1;
//    if (column_flag) wgtflag = 3;//

//    // allocate the memory for the global node numbers
//    int *nnode_pe = (int *) tl_alloc(sizeof(int), npes);
//    idx_t *vtxdist = (idx_t *) tl_alloc(sizeof(idx_t), npes + 1); // the number of nodes belonging to each pe
//    
//    ierr_code = MPI_Allgather(&(local_nnode), 1, MPI_INT, nnode_pe, 1, MPI_INT, g->smpi->ADH_COMM);
//    if (ierr_code != MPI_SUCCESS) messg_err(ierr_code);
//    
//    // gather number of actual nodes on each processor to all processors
//    vtxdist[0] = 0;
//    for (i = 1; i <= npes; i++){
//        vtxdist[i] = vtxdist[i - 1] + nnode_pe[i - 1];
//    }
//    
//    // graph node weights
//    idx_t *vwgt = NULL;
//    if (column_flag) {
//        vwgt = (idx_t *) tl_alloc(sizeof(idx_t), local_nnode);
//        weight_graph_nodes(g, local_nnode, vwgt);
//    }
//    
//    idx_t ncon=1; //number of weights per vertex
//    real_t *tpwgts = (real_t *) tl_alloc(sizeof(real_t), npes);  // fraction of vertex weight to each part ncon x npart array
//    real_t *ubvec  = (real_t *) tl_alloc(sizeof(real_t), ncon);  // imbalance tolerance for each ncon weight
//    ubvec[0] = 1.05;
//    for (i = 0; i < npes; i++) tpwgts[i] = 1.0/npes;//

//    // initialize the hash table
//    //tag(MPI_COMM_WORLD);
//    edge_hashtab_init(edge_hashtab);
//    //tag(MPI_COMM_WORLD);
//    
//    // load nodal connections in hash table to remove redundancies in element connections
//    if (g->ndim == 2) { // 2D
//        edge_hashtab_load_2D(edge_hashtab,g,NULL);
//    } else { // 3D
//        if (column_flag) {
//            edge_hashtab_load_columnar(edge_hashtab,g,NULL);
//        } else {
//            edge_hashtab_load_unstructured(edge_hashtab,g,NULL);
//        }
//    }
//    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
//                                           
//    // count total number of connections and number of connections per node
//    int *nedges = (int *) tl_alloc(sizeof(int), local_nnode);     // the number of nodal connections for each node
//    double dist_max = 0;
//    int nedge_total = 0;
//    dist_max = count_node_edges(g,local_nnode,g->nodeID_3d_to_2d_sur,edge_hashtab,nedges,column_flag);
//    
//    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
//    
//    idx_t *xadj = (idx_t *) tl_alloc(sizeof(idx_t), local_nnode + 1); // the beginning of each nodes list of edges
//    for (i = 0, nedge_total = 0, xadj[0] = 0; i < local_nnode; i++) {
//        nedge_total += nedges[i];
//        xadj[i + 1] = nedge_total;
//    }
//    
//    if (DEBUG == 1) printf("partition_form: pe %d: # edges %d\n", myid, nedge_total);
//    
//    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
//    
//    // load the connections in the adjacency list array
//    idx_t *adjncy = (idx_t *) tl_alloc(sizeof(idx_t), nedge_total); // the nodal connection table
//    idx_t *adjwgt = (idx_t *) tl_alloc(sizeof(idx_t), nedge_total); // edge weights for the graph
//    store_node_edges(g,local_nnode,dist_max,edge_hashtab,nedges,adjncy,adjwgt,xadj,column_flag,NULL);
//    
//    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
//    
//    // partitioning
//    idx_t nparts = npes;
//    idx_t *metis_part = (idx_t *) tl_alloc(sizeof(idx_t), local_nnode); // the partition returned from metis
//    MPI_Comm comm;
//    MPI_Comm_dup(g->smpi->ADH_COMM, &comm);
//    
//    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
//    
//    /* CJT :: debug
//     char fn[30+1];
//     snprintf(fn, 30, "partition_info_pe%d.dat", g->smpi->myid);
//     FILE *fp = fopen(fn, "w");
//     fprintf(fp,"local_nnode: %d \t nedge_total: %d\n",local_nnode,nedge_total);
//     for (i=0; i<nedge_total; i++) {fprintf(fp,"edge: %d \t adjncy: %d \t adjwgt: %d\n",i,adjncy[i],adjwgt[i]);}
//     for (i=0; i<local_nnode; i++) {fprintf(fp,"local_nnode: %d \t vwgt: %d \t xadj: %d\n",i,vwgt[i],xadj[i]);}
//     for (i=0; i<npes+1; i++) {fprintf(fp,"pe: %d \t vtxdist: %d \n",i,vtxdist[i]);}
//     for (i=0; i<npes; i++) {fprintf(fp,"pe: %d \t nnode_pe: %d\n",i,nnode_pe[i]);}
//     fflush(fp);
//     fclose(fp);
//     tl_error("for now");
//     */
//    
//    idx_t edgecut;                       // number of edge cuts (returned)
//    idx_t numflag = 0;                   // flag to indicate numbering scheme
//    idx_t options[4] = { 0, 2, 150, 0 }; // the options for metis
//    ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, &edgecut, metis_part, &comm);
//    
//    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
//    
//    MPI_Comm_free(&comm);
//    
//    
//    if (DEBUG) printf("partition_form: nparts = %d, edgecut = %d\n", nparts, edgecut);
//    
//    // cast the metis partition into part
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
//    
//    metis_part = (idx_t *) tl_free(sizeof(idx_t), local_nnode, metis_part);
//    nnode_pe =     (int *) tl_free(sizeof(int), npes, nnode_pe);
//    vtxdist =    (idx_t *) tl_free(sizeof(idx_t), npes + 1, vtxdist);
//    xadj =       (idx_t *) tl_free(sizeof(idx_t), local_nnode + 1, xadj);
//    nedges =     (idx_t *) tl_free(sizeof(idx_t), local_nnode, nedges);
//    adjncy =     (idx_t *) tl_free(sizeof(idx_t), nedge_total, adjncy);
//    adjwgt =     (idx_t *) tl_free(sizeof(idx_t), nedge_total, adjwgt);
//    tpwgts =    (real_t *) tl_free(sizeof(real_t), npes, tpwgts);
//    ubvec =     (real_t *) tl_free(sizeof(real_t), ncon, ubvec);
//    
//    if (vwgt != NULL) vwgt = (idx_t *) tl_free(sizeof(idx_t), local_nnode, vwgt);
//    tl_list_free_all(EDGE_LIST);
//}