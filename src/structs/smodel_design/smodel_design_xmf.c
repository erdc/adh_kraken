/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  This file contains functions for writing SMODEL\_DESIGN variables to xdmf
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
static int DEBUG_WITH_PICKETS = OFF;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
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
 * \note
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
        
        sol   = dm->superModel[isup].sol;
        ip    = &dm->superModel[isup].ivar_pos;
        ivars = dm->superModel[isup].ivars;
        
        // nodal depth
        if (ip->h != UNSET_INT) {
            data_extract(sol,ivars,ip->h,n,data); // TAKE OUT JUST FOR TESTING!!!!
            if (dm->ts.nt > 0) sarray_scale_replace_dbl(data, dm->ts.nt * 10, 3*n);
            xmf_write_ts_attribute(g, dm->xmf, dm->filebase, ip->name.h,1,dm->ts.nt,true);
            hdf5_write_data(g, dm->filebase, ip->name.h, data, 1, true, dm->ts.nt);
            sarray_init_dbl(data,3*n);
        }

        // nodal velocity
        if (ip->w != UNSET_INT && ip->v != UNSET_INT && ip->u != UNSET_INT) {
            data_extract_vec3d(sol,ivars,ip->u,ip->v,ip->w,n,data);
            xmf_write_ts_attribute(g, dm->xmf, dm->filebase, "velocity",3,dm->ts.nt,true);
            hdf5_write_data(g, dm->filebase, "velocity", data, 3, true, dm->ts.nt);
            sarray_init_dbl(data,3*n);
        } else if (ip->v != UNSET_INT && ip->u != UNSET_INT) {
            data_extract_vec2d(sol,ivars,ip->u,ip->v,n,data);
            xmf_write_ts_attribute(g, dm->xmf, dm->filebase, "velocity",3,dm->ts.nt,true);
            hdf5_write_data(g, dm->filebase, "velocity", data, 3, true, dm->ts.nt);
            sarray_init_dbl(data,3*n);
        } else if (ip->u != UNSET_INT) {
            data_extract(sol,ivars,ip->u,n,data);
            xmf_write_ts_attribute(g, dm->xmf, dm->filebase, "velocity",1,dm->ts.nt,true);
            hdf5_write_data(g, dm->filebase, "velocity", data, 1, true, dm->ts.nt);
            sarray_init_dbl(data,3*n);
        }
        
        // nodal depth-averaged velocity
        if (ip->vda != UNSET_INT && ip->uda != UNSET_INT) {
            data_extract_vec2d(sol,ivars,ip->uda,ip->vda,n,data);
            xmf_write_ts_attribute(g, dm->xmf, dm->filebase, "da-velocity",3,dm->ts.nt,true);
            hdf5_write_data(g, dm->filebase, "da-velocity", data, 3, true, dm->ts.nt);
            sarray_init_dbl(data,3*n);
        } else if (ip->uda != UNSET_INT) {
            data_extract(sol,ivars,ip->uda,n,data);
            xmf_write_ts_attribute(g, dm->xmf, dm->filebase, "da-velocity",1,dm->ts.nt,true);
            hdf5_write_data(g, dm->filebase, "da-velocity", data, 1, true, dm->ts.nt);
            sarray_init_dbl(data,3*n);
        }
        
        // constituents
        for (int icon=0; icon<ip->ntrns; icon++) {
            if (ip->con[icon] != UNSET_INT) {
                data_extract(sol,ivars,ip->con[icon],n,data);
                xmf_write_ts_attribute(g, dm->xmf, dm->filebase, ip->name.con[icon],1,dm->ts.nt,true);
                hdf5_write_data(g, dm->filebase, ip->name.con[icon], data, 1, true, dm->ts.nt);
                sarray_init_dbl(data,3*n);
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
