/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  This fil initializes and SMODEL_DESIGN structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = ON;
static int DEBUG_WITH_PICKETS = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes an AdH Designer Model
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] dmod           (SDMODEL **)  a pointer to an AdH design-level model
 * @param[in] argc            (int *) the number of command line arguements
 * @param[in] argv            (int *) the command line arguements

 * \note CJT :: still need to (1) Prep PETSc, (2) Open output and write initial Conditions, etc.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int smodel_design_init(SMODEL_DESIGN *dmod, char *filebase, bool input_check
#ifdef _MESSG
                   ,MPI_Comm *comm
#endif
) {
    FILE *fp = NULL;
    int ierr = 0, imono = 0;
    
    // Is the file supplied a super file?
    char fname[MAXLINE];
    strcpy(fname,filebase);
    strcat(fname,".design");
    if (DEBUG) {
        printf("------------------------------------------------------\n");
        printf("------------------------------------------------------\n");
        printf("Initiating Designer Model '%s' Read\n",fname);
        printf("------------------------------------------------------\n");
    }
    if (sfile_checkExists(fname)) {
        smodel_design_read(dmod,fname);
    } else {
       //printf("single file read\n");
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Read the grid
    //++++++++++++++++++++++++++++++++++++++++++++++
    if (DEBUG) {
        printf("------------------------------------------------------\n");
        printf("------------------------------------------------------\n");
        printf("Initiating Designer Model Grid Read\n");
        printf("------------------------------------------------------\n");
    }
#ifdef _MESSG
    sgrid_read(&(dmod->grid),dmod->filename_grid,model_comm);
#else
    sgrid_read(&(dmod->grid),dmod->filename_grid);
#endif
    if (DEBUG) {
        sgrid_printScreen(dmod->grid);
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Assign designer grid to each superModel
    //++++++++++++++++++++++++++++++++++++++++++++++
    for (imono=0; imono<dmod->nSuperModels; imono++) {
        dmod->superModel[imono].grid = dmod->grid;
        //smodel_super_read(&(dmod->superModel[imono]))
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Read superModel parameters
    //++++++++++++++++++++++++++++++++++++++++++++++
    for (imono=0; imono<dmod->nSuperModels; imono++) {
        if (DEBUG) {
            printf("------------------------------------------------------\n");
            printf("------------------------------------------------------\n");
            printf("Initiating Super Model '%s' File Read\n",dmod->superModel[imono].filebase);
            printf("------------------------------------------------------\n");
        }
        smodel_super_read(&(dmod->superModel[imono]));
        //exit(-1);

        if (DEBUG) {
            printf("------------------------------------------------------\n");
            printf("Reading SuperModel '%s' Initialization File\n",dmod->superModel[imono].filebase);
            printf("------------------------------------------------------\n");
        }
        smodel_super_read_init(&(dmod->superModel[imono]),dmod->superModel[imono].filebase);
    }

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Open SuperModel Parameter Coverage File
    // One coverage can be used for multiple superModels
    //++++++++++++++++++++++++++++++++++++++++++++++filebase,NULL,NULL,suffix,mode,TRUE);
    if (DEBUG) {
        printf("------------------------------------------------------\n");
        printf("------------------------------------------------------\n");
        printf("Initiating Designer Model Parameter File Read\n");
        printf("------------------------------------------------------\n");
    }
    scoverage_read(&dmod->params,filebase);
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Allocate linear systems
    //++++++++++++++++++++++++++++++++++++++++++++++



    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Stop code if only file checking
    if (input_check) {tl_error("Exiting after file check!");}
    //++++++++++++++++++++++++++++++++++++++++++++++




    return ierr;
}

