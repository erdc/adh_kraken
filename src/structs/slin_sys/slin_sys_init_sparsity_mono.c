#include "adh.h"
static int **temp_cols_diag = NULL;
static int **temp_cols_off_diag = NULL;
static int *nnz_rows_diag = NULL;
static int *nnz_rows_off_diag = NULL;
static int *nnz_rows_diag_no_duplicate = NULL;
static int *nnz_rows_off_diag_no_duplicate = NULL;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function allocates the sparsity pattern of the split CSR structure,
 *  it allocates and sets the indptr, and cols of diag and off diag given the physics material
 *  layouts. It also determines number of nonzeros by allocating the vals.
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *   @param[in,out] sm (SLIN_SYS*) - pointer to an instant of the SLIN_SYS struct,
 *  this will contain pointer to the CSR structure
 *  @param[in] sm (SMODEL_SUPER*) - pointer to an instant of the SMODEL_SUPER struct,
 *  this will contain pointer to grid, physics materials
 * \note For now, we are doing the most memory safe but computationally redundant way. 
 * There is some potential to speed up the routine by guessing number of nonzeros per row.
 * For now there is no guessing.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void slin_sys_init_sparsity_mono(SLIN_SYS *lin_sys, int *elem3d_physics_mat_id, 
    int *elem2d_physics_mat_id, int *elem1d_physics_mat_id, SMAT_PHYSICS *elem_physics_mat,
    SGRID *grid, int **ivars){

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Creates sparsity format in split CSR format
    // Note: this amounts to establishing graph of the degrees of freedom of our system
    //
    //now with info from the super models, assign appropriate vals to the design model
    //for each unique model we build up the sparsity and solution variables
    //allocate CSR sparsity structure for each unique model  
    //missing a step to find proper index of supermodel that is unique
    //THIS SHOULDNT BE i, this should be the ith non-simple super model
    //temporarily stores column positions of diagonal and off diagonal blocks
    //diagonal blocks are local column numbers, off-diagonal are global
    //++++++++++++++++++++++++++++++++++++++++++++++

    bool has_off_diag = false;
    int row,col;
    int nnodes;
    int nvars_elem;
    int ndofs_ele;
    int i,j,k,l, mat_id;
    int elem_vars[adh_def.n_ivars];
    int dofs[adh_def.n_ivars*MAX_NNODE],global_dofs[adh_def.n_ivars*MAX_NNODE];
    int isize_prev;
    int NNZ_diag = 0;
    int nnz_row_diag;
    int NNZ_off_diag=0;
    int nnz_row_off_diag;
    int count_diag=0;
    int *rn_diag;
    int count_off_diag = 0;
    int *rn_off_diag;
    int *ghosts = lin_sys->ghosts;

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // We know that each row contains at least one non-zero entry, 
    // and the number of rows is known upfront. 
    // Thus, (i,j)-pairs can be represented as a vector of vectors
    //first allocate this by findinf max size of nnz per row
    //++++++++++++++++++++++++++++++++++++++++++++++


    int *local_range = lin_sys->local_range;
    int nrows = *(lin_sys->local_size);
    int isize = *(lin_sys->local_size_old);
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // see if we have off diagonal blocks to handle
    //++++++++++++++++++++++++++++++++++++++++++++++

    if(lin_sys->indptr_off_diag!=NULL){
        has_off_diag=true;
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // reallocate static variables if necessary
    //++++++++++++++++++++++++++++++++++++++++++++++

    if (isize < nrows){
        isize_prev = isize;
        isize = nrows;

        nnz_rows_diag = (int *) tl_realloc(sizeof(int), isize, isize_prev, nnz_rows_diag);
        nnz_rows_diag_no_duplicate = (int *) tl_realloc(sizeof(int), isize, isize_prev, nnz_rows_diag_no_duplicate);
        if (has_off_diag){
            nnz_rows_off_diag = (int *) tl_realloc(sizeof(int), isize, isize_prev, nnz_rows_off_diag);
            nnz_rows_off_diag_no_duplicate = (int *) tl_realloc(sizeof(int), isize, isize_prev, nnz_rows_off_diag_no_duplicate);
        }
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // update variables in lin sys object
    //++++++++++++++++++++++++++++++++++++++++++++++
    *(lin_sys->local_size_old) = isize;

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // initialize arrays with zeros
    //++++++++++++++++++++++++++++++++++++++++++++++

    sarray_init_int(nnz_rows_diag, nrows);
    if (has_off_diag){
        sarray_init_int(nnz_rows_off_diag, nrows);
    }

    
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // First set of loops is solely to establish how many nonzeros are there 
    // and how to dynamically allocate temporary sparsity arrays
    //++++++++++++++++++++++++++++++++++++++++++++++

    for (j=0;j<grid->nelems3d;j++){

        nnodes = grid->elem3d[j].nnodes;
        //pull all global information to local memory
        mat_id = elem3d_physics_mat_id[j];

        // Get stuff from physics mat
        nvars_elem = elem_physics_mat[mat_id].ivar_pos.n;

        ndofs_ele = nnodes*nvars_elem;

        //++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++
        // for residual we only need dof numbers local to process (including ghost nodes)
        // this is a complicated map but maybe we can simplify in simpler cases by 
        // replacing different routine
        // usually would take the local cell number and compute the associated dofs
        // but this has expanded arguments so it will work for elem1d,elem2d as well
        // cell # is implicit
        //++++++++++++++++++++++++++++++++++++++++++++++
        get_cell_dofs_ivars(dofs, ivars, nnodes, grid->elem3d[j].nodes, nvars_elem, elem_physics_mat[mat_id].ivar_loc);
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
    
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (k=0;k<ndofs_ele;k++){
            //current row number local to the process
            row = dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
            if(row<nrows){
                for (l=0;l<ndofs_ele;l++){
                    //differentiate  between diag and off diag blocks by column number (local to process)
                    col = dofs[l];
                    if (col<nrows){
                        //then this is on diagonal block and we want column numbers local to proces
                        nnz_rows_diag[row]+=1;
                    }
                    else{
                        //otherwise we are on off-diagonal and need global dof
                        nnz_rows_off_diag[row]+=1;//
                    }
                }
            }
        }
    }
    
    for (j=0;j<grid->nelems2d;j++){
        nnodes = grid->elem2d[j].nnodes;
        //pull all global information to local memory
        mat_id = elem2d_physics_mat_id[j];

        // Get stuff from physics mat
        nvars_elem = elem_physics_mat[mat_id].ivar_pos.n;

        ndofs_ele = nnodes*nvars_elem;

        //++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++
        // for residual we only need dof numbers local to process (including ghost nodes)
        // this is a complicated map but maybe we can simplify in simpler cases by 
        // replacing different routine
        // usually would take the local cell number and compute the associated dofs
        // but this has expanded arguments so it will work for elem1d,elem2d as well
        // cell # is implicit
        //++++++++++++++++++++++++++++++++++++++++++++++
        get_cell_dofs_ivars(dofs, ivars, nnodes, grid->elem2d[j].nodes, nvars_elem, elem_physics_mat[mat_id].ivar_loc);
        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
    
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (k=0;k<ndofs_ele;k++){
            //current row number local to the process
            row = dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
            if(row<nrows){
                for (l=0;l<ndofs_ele;l++){
                    //differentiate  between diag and off diag blocks by column number (local to process)
                    col = dofs[l];
                    if (col<nrows){
                        //then this is on diagonal block and we want column numbers local to proces
                        nnz_rows_diag[row]+=1;
                    }
                    else{
                        //otherwise we are on off-diagonal and need global dof
                        nnz_rows_off_diag[row]+=1;//
                    }
                }
            }
        }
    }

    for (j=0;j<grid->nelems1d;j++){

        nnodes = grid->elem1d[j].nnodes;
        //pull all global information to local memory
        mat_id = elem1d_physics_mat_id[j];

        // Get stuff from physics mat
        nvars_elem = elem_physics_mat[mat_id].ivar_pos.n;

        ndofs_ele = nnodes*nvars_elem;

        //++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++
        // for residual we only need dof numbers local to process (including ghost nodes)
        // this is a complicated map but maybe we can simplify in simpler cases by 
        // replacing different routine
        // usually would take the local cell number and compute the associated dofs
        // but this has expanded arguments so it will work for elem1d,elem2d as well
        // cell # is implicit
        //++++++++++++++++++++++++++++++++++++++++++++++
        get_cell_dofs_ivars(dofs, ivars, nnodes, grid->elem1d[j].nodes, nvars_elem, elem_physics_mat[mat_id].ivar_loc);

        //this gets global dofs from local dofs, and fmapglobal is this best way to do it?
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
    
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (k=0;k<ndofs_ele;k++){
            //current row number local to the process
            row = dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
            if(row<nrows){
                for (l=0;l<ndofs_ele;l++){
                    //differentiate  between diag and off diag blocks by column number (local to process)
                    col = dofs[l];
                    if (col<nrows){
                        //then this is on diagonal block and we want column numbers local to proces
                        nnz_rows_diag[row]+=1;
                    }
                    else{
                        //otherwise we are on off-diagonal and need global dof
                        nnz_rows_off_diag[row]+=1;//
                    }
                }
            }
        }
    }


    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Use nnz_rows to dynamically allocate
    //++++++++++++++++++++++++++++++++++++++++++++++
    
   temp_cols_diag = (int**) tl_alloc(sizeof(int*), nrows);
    for(j=0;j<nrows;j++){
        temp_cols_diag[j] = (int*) tl_alloc(sizeof(int), nnz_rows_diag[j]);
        for(k=0;k<nnz_rows_diag[j];k++){
            temp_cols_diag[j][k]=INT_MAX;
        }
    }
    if (has_off_diag){
            temp_cols_off_diag = (int**) tl_alloc(sizeof(int*), nrows);
        for(j=0;j<nrows;j++){
            temp_cols_off_diag[j] = (int*) tl_alloc(sizeof(int), nnz_rows_off_diag[j]);
            for(k=0;k<nnz_rows_off_diag[j];k++){
                temp_cols_off_diag[j][k]=INT_MAX;
            }
        }   
    }

    //Seems redundant but must reuse as indexing
    sarray_init_int(nnz_rows_diag, nrows);
    if (has_off_diag){
        sarray_init_int(nnz_rows_off_diag, nrows);
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Loop thorugh each element and find sparsity
    //++++++++++++++++++++++++++++++++++++++++++++++
    
    for (j=0;j<grid->nelems3d;j++){
        nnodes = grid->elem3d[j].nnodes;
        mat_id = elem3d_physics_mat_id[j];
        nvars_elem = elem_physics_mat[mat_id].ivar_pos.n;
        ndofs_ele = nnodes*nvars_elem;
        get_cell_dofs_ivars(dofs, ivars, nnodes, grid->elem3d[j].nodes, nvars_elem, elem_physics_mat[mat_id].ivar_loc);
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
        
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (k=0;k<ndofs_ele;k++){
            //current row number local to the process
            row = dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
            if(row<nrows){
                for (l=0;l<ndofs_ele;l++){
                    //differentiate  between diag and off diag blocks by column number (local to process)
                    col = dofs[l];
                    if (col<nrows){
                        //then this is on diagonal block and we want column numbers local to proces
                        temp_cols_diag[row][nnz_rows_diag[row]]=col;
                        nnz_rows_diag[row]+=1;
                    }
                    else{
                        //otherwise we are on off-diagonal and need global dof
                        temp_cols_off_diag[row][nnz_rows_off_diag[row]]=global_dofs[l];
                        nnz_rows_off_diag[row]+=1;
                    }
                }
            }
        }
    }

    for (j=0;j<grid->nelems2d;j++){
        nnodes = grid->elem2d[j].nnodes;
        mat_id = elem2d_physics_mat_id[j];
        nvars_elem = elem_physics_mat[mat_id].ivar_pos.n;
        ndofs_ele = nnodes*nvars_elem;
        get_cell_dofs_ivars(dofs, ivars, nnodes, grid->elem2d[j].nodes, nvars_elem, elem_physics_mat[mat_id].ivar_loc);
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
        
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (k=0;k<ndofs_ele;k++){
            //current row number local to the process
            row = dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
            if(row<nrows){
                for (l=0;l<ndofs_ele;l++){
                    //differentiate  between diag and off diag blocks by column number (local to process)
                    col = dofs[l];
                    if (col<nrows){
                        //then this is on diagonal block and we want column numbers local to proces
                        temp_cols_diag[row][nnz_rows_diag[row]]=col;
                        nnz_rows_diag[row]+=1;
                    }
                    else{
                        //otherwise we are on off-diagonal and need global dof
                        temp_cols_off_diag[row][nnz_rows_off_diag[row]]=global_dofs[l];
                        nnz_rows_off_diag[row]+=1;
                    }
                }
            }
        }
    }

    for (j=0;j<grid->nelems1d;j++){
        nnodes = grid->elem1d[j].nnodes;
        mat_id = elem1d_physics_mat_id[j];
        nvars_elem = elem_physics_mat[mat_id].ivar_pos.n;
        ndofs_ele = nnodes*nvars_elem;
        get_cell_dofs_ivars(dofs, ivars, nnodes, grid->elem1d[j].nodes, nvars_elem, elem_physics_mat[mat_id].ivar_loc);
        local_dofs_to_global_dofs(global_dofs,ndofs_ele,dofs,local_range,ghosts);
        
        //keep sparsity pattern in temp_cols_diag, temp_cols_off_diag
        for (k=0;k<ndofs_ele;k++){
            //current row number local to the process
            row = dofs[k];
            //loop through each row of the local elemental matrix
            //we only need to fill in columns if the current row is local to process
            if(row<nrows){
                for (l=0;l<ndofs_ele;l++){
                    //differentiate  between diag and off diag blocks by column number (local to process)
                    col = dofs[l];
                    if (col<nrows){
                        //then this is on diagonal block and we want column numbers local to proces
                        temp_cols_diag[row][nnz_rows_diag[row]]=col;
                        nnz_rows_diag[row]+=1;
                    }
                    else{
                        //otherwise we are on off-diagonal and need global dof
                        temp_cols_off_diag[row][nnz_rows_off_diag[row]]=global_dofs[l];
                        nnz_rows_off_diag[row]+=1;
                    }
                }
            }
        }
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // now that all nnz_cols* and nnz* are filled, we need to sort and remove duplicates
    // counter for the total number of non-zero entries
    // note nnz-rows is not actually nnz per row, includes duplicates
    // could be done with hybrid omp too?
    // do for diagonal and off diagonal blocks
    //++++++++++++++++++++++++++++++++++++++++++++++

    for (i=0;i<nrows;i++){
        // sort the column indices (j-entries)
        //use stdlib.h qsort
        qsort(temp_cols_diag[i], nnz_rows_diag[i], sizeof(int), compare_ints);
        //this should remove duplicates
        nnz_row_diag = sarray_unique_int(temp_cols_diag[i], nnz_rows_diag[i]);
        //overwrite nnz row with sarray_unique_int
        nnz_rows_diag_no_duplicate[i] = nnz_row_diag;
        //add nnz in a row to the NNZ
        //maybe overwrire nnz_rows[i] instead if we want this stored, and then sum it
        NNZ_diag+=nnz_row_diag;
    }

    lin_sys->nnz_diag = NNZ_diag;

    // off diagonalblocks
    if (has_off_diag){
        for (i=0;i<nrows;i++){
            // sort the column indices (j-entries)
            //use stdlib.h qsort
            qsort(temp_cols_off_diag[i], nnz_rows_off_diag[i], sizeof(int), compare_ints);
            //this should hopefully remove duplicates?
            nnz_row_off_diag = sarray_unique_int(temp_cols_off_diag[i], nnz_rows_off_diag[i]);
            //overwrite nnz row with sarray_unique_int?
            nnz_rows_off_diag_no_duplicate[i] = nnz_row_off_diag;
            //add nnz in a row to the NNZ
            //maybe overwrire nnz_rows[i] instead if we want this stored, and then sum it
            NNZ_off_diag+=nnz_row_off_diag;
        }
        lin_sys->nnz_off_diag = NNZ_off_diag;
    }

    // resize as needed
    lin_sys->cols_diag = (int *) tl_realloc(sizeof(int), lin_sys->nnz_diag, lin_sys->nnz_diag_old, lin_sys->cols_diag);
    lin_sys->vals_diag = (double *) tl_realloc(sizeof(double), lin_sys->nnz_diag, lin_sys->nnz_diag_old, lin_sys->vals_diag);
    if (has_off_diag){
        lin_sys->cols_off_diag = (int *) tl_realloc(sizeof(int), lin_sys->nnz_off_diag, lin_sys->nnz_off_diag_old, lin_sys->cols_off_diag);
        lin_sys->vals_off_diag = (double *) tl_realloc(sizeof(double), lin_sys->nnz_off_diag, lin_sys->nnz_off_diag_old, lin_sys->vals_off_diag);
    }
    //now use info to fill in indptr and cols
    for(i=0;i<nrows;i++){
        lin_sys->indptr_diag[i] = count_diag;
        rn_diag = temp_cols_diag[i];
        //think about how to do this since each temp_rows may be different size
        for(j=0;j<nnz_rows_diag_no_duplicate[i];j++){
            lin_sys->cols_diag[count_diag] = rn_diag[j];
            count_diag++;
        }
    }
    //also last sm->indptr  needs to be the last entry
    lin_sys->indptr_diag[nrows] = count_diag;    
    if (has_off_diag){
        for(i=0;i<nrows;i++){
            lin_sys->indptr_off_diag[i] = count_off_diag;
            rn_off_diag = temp_cols_off_diag[i];
            //think about how to do this since each temp_rows may be different size
            for(j=0;j<nnz_rows_off_diag_no_duplicate[i];j++){
                lin_sys->cols_off_diag[count_off_diag] = rn_off_diag[j];
                count_off_diag++;
            }
        }
        lin_sys->indptr_off_diag[nrows] = count_off_diag;
    }



    //for now free calls are in script, can separate to new routine later
    for (i=0; i<nrows; i++) {
        temp_cols_diag[i] = (int*) tl_free(sizeof(int), nnz_rows_diag[i], temp_cols_diag[i]);
    }
    temp_cols_diag = (int **) tl_free(sizeof(int *), nrows, temp_cols_diag);
    nnz_rows_diag= (int *) tl_free(sizeof(int), nrows, nnz_rows_diag);
    nnz_rows_diag_no_duplicate= (int *) tl_free(sizeof(int), nrows, nnz_rows_diag_no_duplicate);

    if (has_off_diag){
        for (i=0; i<nrows; i++) {
            temp_cols_off_diag[i] = (int*) tl_free(sizeof(int), nnz_rows_off_diag[i], temp_cols_off_diag[i]);
        }
        temp_cols_off_diag = (int **) tl_free(sizeof(int *), nrows, temp_cols_off_diag);
        nnz_rows_off_diag= (int *) tl_free(sizeof(int), nrows, nnz_rows_off_diag);
        nnz_rows_off_diag_no_duplicate= (int *) tl_free(sizeof(int), nrows, nnz_rows_off_diag_no_duplicate);
    }
    
    printf("CSR sparsity completed\n");
}
