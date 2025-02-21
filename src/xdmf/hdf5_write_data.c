
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This file write data to an HDF5 file
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] g                (SGRID *)  an AdH grid
 * @param[in] fbase       (char *)  the hdf5 filename base
 * @param[in] name         (char *)  the variable name
 * @param[in] data         (double *) the data to be written
 * @param[in] ndim         (int)   the variable dimensionality
 * @param[in] isNodal  (bool)  is the data on nodes - true or elements - false
 * @param[in] nt             (int) the AdH time-step
 
 *
 * \note cjt - only works for serial at the moment
 * \note Tutorial at: https://docs.hdfgroup.org/hdf5/v1_14/v1_14_4/_ex_a_p_i.html
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#define RANK 2
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void hdf5_write_data(SGRID *g, char *fbase, char *name, double *data, int ndim, bool isNodal, int nt) {
    
    //++++++++++++++++++++++++++++++++++++++
    // Write HD5 Data
    //++++++++++++++++++++++++++++++++++++++
    hid_t file_id, dset_id, memspace, filespace, plist_id, grp;
    hsize_t dims[RANK], count[RANK], offset[RANK];
    herr_t status;
    
    char number[50] = ""; sprintf(number, "%d", nt);
    
    //++++++++++++++++++++++++++++++++++++++
    //Open file for parallel i/o
    //++++++++++++++++++++++++++++++++++++++
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
#ifdef _MESSG
    MPI_Info info  = MPI_INFO_NULL;
    H5Pset_fapl_mpio(plist_id, g->smpi->ADH_COMM, info);
#endif
    char fname[50];
    strcpy(fname,fbase); strcat(fname, ".h5");
    file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);
    
    //++++++++++++++++++++++++++++++++++++++
    // Either create (firstWrite) or open group
    //++++++++++++++++++++++++++++++++++++++
    char path[50] = "/Data/";
    strcat(path,name);
    if (nt == 0) {
        grp = H5Gcreate(file_id, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    } else {
        grp = H5Gopen(file_id, path, H5P_DEFAULT);
    }
        
    if (isNodal) {
        //++++++++++++++++++++++++++++++++++++++
        // Node-based data
        //++++++++++++++++++++++++++++++++++++++
        dims[0] = g->macro_nnodes;
        count[0] = g->my_nnodes;
        offset[0] = 0; // mpi_rank * count[0]; // cjt -- serial only
    } else {
        //++++++++++++++++++++++++++++++++++++++
        // Element-based data
        //++++++++++++++++++++++++++++++++++++++
        dims[0] = g->macro_nSegs +  g->macro_nQuads + g->macro_nTris + g->macro_nTets + g->macro_nPrisms;
        count[0] = g->my_nSegs + g->my_nQuads + g->my_nTris + g->my_nTets + g->my_nPrisms;
        offset[0] = 0; // mpi_rank * count[0]; // cjt -- serial only
    }
    
    dims[1]   = ndim;
    count[1]  = ndim;  // # of data columns = data dimension
    offset[1] = 0;     // for hpc later
    
    //sarray_printScreen_dbl(data,dims[0]*dims[1],"sol");
    
    //filespace = H5Screate_simple(dims[1], dims, NULL);
    filespace = H5Screate_simple(RANK, dims, NULL);
    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp, number, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    // give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    //memspace  = H5Screate_simple(dims[1], count, NULL);
    memspace  = H5Screate_simple(RANK, count, NULL);
    
    //determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    //Create property list for collective dataset write.
    //declare collextive data file weiting
    plist_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef _MESSG
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif
    
    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);


    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp);
    
    //++++++++++++++++++++++++++++++++++++++
    H5Fclose(file_id);
    //++++++++++++++++++++++++++++++++++++++
    
    if (DEBUG) {
        if (g->smpi->myid == 0) {printf("PHDF5 example finished with no errors\n");}
    }
}
