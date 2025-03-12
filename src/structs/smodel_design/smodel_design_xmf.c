/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  This file contains functions for writing SMODEL\_DESIGN variables to xdmf
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
static int DEBUG_WITH_PICKETS = OFF;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void data_stack(double **matrix, int ivar, int ndim, int n, double *data);
void data_extract(double *sol, int **ivars, int varpos, int n, double *data);
void data_extract_vec2d(double *sol, int **ivars, int varpos1, int varpos2, int n, double *data);
void data_extract_vec3d(double *sol, int **ivars, int varpos1, int varpos2, int varpos3, int n, double *data);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Initializes a XMF/HDF5 output and writes the initial grid
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] dm  (SMODEL_DESIGN *)  a pointer to a design model structure
 * @param[in] filename  (char *)  the base xmf/h5 filename
 * @param[in] domain_name  (char *)  the xmf AdH domain name
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_design_xmf_init(SMODEL_DESIGN *dm, char *filename, char *domain_name) {
    dm->xmf = xmf_init(filename,domain_name,filename);
    assert(dm->xmf != NULL);
    hdf5_init(dm->grid, filename);       // initialize HDF5
    hdf5_write_grid(dm->grid, filename); // write grid file
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Print SMODEL\_DESIGN data in HDF5/XDMF format
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] dm  (SMODEL_DESIGN *)  a pointer to a design model structure
 * @param[in] mesh_no (int) the adapted mesh ID.  If not printing adapted meshes, this is always 0;
 *
 * \note CJT -- figure out later how to write just u for 1D SW for example
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_design_xmf_write(SMODEL_DESIGN *dm, int mesh_no) {
#ifdef _HDF5
    int i,j;
    
    //++++++++++++++++++++++++++++++
    // aliases
    //++++++++++++++++++++++++++++++
    SGRID *g = dm->grid;
    int myid = g->smpi->myid;
    double *sol = NULL;
    SIVAR_POSITION *ip = NULL;
    int **ivars = NULL;
    int n = g->my_nnodes;
    SMODEL_SUPER *sm;
    
    //++++++++++++++++++++++++++++++
    // data buffer
    //++++++++++++++++++++++++++++++
    double data[3*n]; sarray_init_dbl(data,3*n);
    
    //++++++++++++++++++++++++++++++
    // prep grid name for data
    //++++++++++++++++++++++++++++++
    char gridname[50] = "";
    strcpy(gridname,dm->filebase);
    char number[50] = "";
    sprintf(number, "%d", dm->ts.nt);
    strcat(gridname,number);
    
    //++++++++++++++++++++++++++++++
    // write xmf initial time-step matter
    //++++++++++++++++++++++++++++++
    xmf_write_ts_header(myid,dm->xmf,gridname,dm->ts.time);
    xmf_write_geometry(g,dm->xmf,dm->filebase,mesh_no);
    xmf_write_topology(g,dm->xmf,dm->filebase,mesh_no);
    
    //++++++++++++++++++++++++++++++
    // write data xmf/hdf5
    //++++++++++++++++++++++++++++++
    for (int isup=0; isup<dm->nSuperModels; isup++) {
        
        // write independent nodal variables
        sm    = &dm->superModel[isup];
        sol   = sm->sol;
        ip    = &(sm->ivar_pos);
        ivars = sm->ivars;
        
        for (int ivar=0; ivar<N_IVARS_TOTAL; ivar++) {
            if (ip->var[ivar] != UNSET_INT) {
                printf(">>> writing XDMF: var_name: %s || var_ndim: %d\n",IVAR_NAME[ivar],IVAR_DIM[ivar]);
                if (IVAR_DIM[ivar] == 1) {
                    data_extract(sol,ivars,ip->var[ivar],n,data);
                    xmf_write_ts_attribute(g, dm->xmf, dm->filebase, IVAR_NAME[ivar], 1, dm->ts.nt, true);
                    hdf5_write_data(g, dm->filebase, IVAR_NAME[ivar], data, 1, true, dm->ts.nt);
                    sarray_init_dbl(data,3*n);
                } else if (IVAR_DIM[ivar] == 2) {
                    data_extract_vec2d(sol,ivars,ip->var[ivar],ip->var[ivar+1],n,data);
                    xmf_write_ts_attribute(g, dm->xmf, dm->filebase, IVAR_NAME[ivar], 2, dm->ts.nt, true);
                    hdf5_write_data(g, dm->filebase, IVAR_NAME[ivar], data, 3, true, dm->ts.nt);
                    sarray_init_dbl(data,3*n);
                    ivar += 1;
                } else if (IVAR_DIM[ivar] == 3) {
                    data_extract_vec3d(sol,ivars,ip->var[ivar],ip->var[ivar+1],ip->var[ivar+2],n,data);
                    xmf_write_ts_attribute(g, dm->xmf, dm->filebase, IVAR_NAME[ivar], 3, dm->ts.nt, true);
                    hdf5_write_data(g, dm->filebase, IVAR_NAME[ivar], data, 3, true, dm->ts.nt);
                    sarray_init_dbl(data,3*n);
                    ivar += 2;
                }
            }
        }
        
        // write dependent nodal variables
        SDVAR_POSITION *dp = &sm->dvars.sdvar_pos_node;
        if (sm->dvars.n_dvar > 0) { 
            for (int ivar=0; ivar<N_DVARS; ivar++) {
                // cjt - need to add print flags here to let user choose what dep-variables are written
                if (dp->var[ivar] != UNSET_INT) {
                    printf(">>> writing XDMF: var_name: %s || var_ndim: %d\n",DVAR_NAME[ivar],DVAR_DIM[ivar]);
                    if (DVAR_DIM[ivar] == 1) {
                        data_stack(sm->dvars.nodal_dvar,dp->var[ivar],1,n,data);
                        xmf_write_ts_attribute(g, dm->xmf, dm->filebase, DVAR_NAME[ivar],1,dm->ts.nt,true);
                        hdf5_write_data(g, dm->filebase, DVAR_NAME[ivar], data, 1, true, dm->ts.nt);
                        sarray_init_dbl(data,3*n);
                    } else if (DVAR_DIM[ivar] == 2) {
                        data_stack(sm->dvars.nodal_dvar,dp->var[ivar],2,n,data);
                        xmf_write_ts_attribute(g, dm->xmf, dm->filebase, DVAR_NAME[ivar],2,dm->ts.nt,true);
                        hdf5_write_data(g, dm->filebase, DVAR_NAME[ivar], data, 3, true, dm->ts.nt);
                        sarray_init_dbl(data,3*n);
                        ivar += 1;
                    } else if (DVAR_DIM[ivar] == 3) {
                        data_stack(sm->dvars.nodal_dvar,dp->var[ivar],3,n,data);
                        xmf_write_ts_attribute(g, dm->xmf, dm->filebase, DVAR_NAME[ivar],3,dm->ts.nt,true);
                        hdf5_write_data(g, dm->filebase, DVAR_NAME[ivar], data, 3, true, dm->ts.nt);
                        sarray_init_dbl(data,3*n);
                        ivar += 2;
                    }
                }
            }
        }
        
    }
    
    //++++++++++++++++++++++++++++++
    // close xmf time-step info
    //++++++++++++++++++++++++++++++
    xmf_write_ts_tail(dm->xmf);
    
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void data_extract(double *sol, int **ivars, int varpos, int n, double *data) {
    for (int i=0; i<n; i++) {
        data[i] = sol[ivars[varpos][i]];
    }
}
void data_extract_vec2d(double *sol, int **ivars, int varpos1, int varpos2, int n, double *data) {
    int k=0;
    for (int i=0; i<n; i++) {
        data[k] = sol[ivars[varpos1][i]]; k++;
        data[k] = sol[ivars[varpos2][i]]; k++;
        data[k] = 0.0; k++;
    }
}
void data_extract_vec3d(double *sol, int **ivars, int varpos1, int varpos2, int varpos3, int n, double *data) {
    int k=0;
    for (int i=0; i<n; i++) {
        data[k] = sol[ivars[varpos1][i]]; k++;
        data[k] = sol[ivars[varpos2][i]]; k++;
        data[k] = sol[ivars[varpos3][i]]; k++;
    }
}

void data_stack(double **matrix, int ivar, int ndim, int n, double *data) {
    int idim,i,j,k=0;
    for (i=0; i<n; i++) {
        for (idim=0; idim<ndim; idim++) {
            data[k] = matrix[ivar + idim][i]; k++;
        }
        if (ndim == 2) { // add a 0 to the third dimension
            data[k] = 0.0; k++;
        }
    }
    
//    k=0;
//    for (i=0; i<n; i++) {
//        for (j=0; j<3; j++) {
//            printf("data[%d]: %f\n",k,data[k]);
//            k++;
//        }
//    }
//    exit(-1);
}

