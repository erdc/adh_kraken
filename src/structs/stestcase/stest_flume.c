/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  stestcase.c Collects routine for the AdH flume test case for testing NB boundary conditions */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief Initializes a flume testcase
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] sm (SMODEL_SUPER *)  a pointer to a superModel structure that stores AdH model results
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void testcase_sw2_flume_init(SMODEL_SUPER *sm) {
    // this is initialized through file read
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief Finalizes a flume testcase
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] sm (SMODEL_SUPER *)  a pointer to a superModel structure that stores AdH model results
 * \note
 */
void testcase_sw2_flume_final(SMODEL_SUPER *sm) {
    // nothing needed
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief Writes to file the flume testcase errors
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] sm (SMODEL_SUPER *)  a pointer to a superModel structure that stores AdH model results
 * \note
 */
void testcase_sw2_flume_write(SMODEL_SUPER *sm) {
    
    double convert_to_rads = 3.141592653589793 / 180.;
    
    // the inflow is 0.1 m/s at a 45 degree angle
    double analytic_vx = 0.1/sqrt(2.);
    double analytic_vy = analytic_vx;
    
    double error_vx = 0., error_vx_max;
    double error_vy = 0., error_vy_max;
    double error_h = 0., error_h_max;
    for (int inode=0; inode<sm->grid->my_nnodes; inode++) {
        error_vx += fabs(get_node_uda(sm,inode) - analytic_vx);
        error_vy += fabs(get_node_vda(sm,inode) - analytic_vy);
        error_h  += fabs(get_node_h(sm,inode) - 1.);  // this is just to see how far from rigid lid
    }
    
//    double grid_mass_error = tl_find_grid_mass_error_elem2d(mod->density, mod->sw->d2->head, mod->sw->d2->vel, sm->grid, mod->flag, mod->initial_grid_mass, mod->series_head, mod->str_values, mod->dt, &total_time_mass_flux);
    
#ifdef _MESSG
    error_vx_max=messg_dsum(error_vx,sm->grid->smpi->ADH_COMM);
    error_vy_max=messg_dsum(error_vy,sm->grid->smpi->ADH_COMM);
    error_h_max=messg_dsum(error_h,sm->grid->smpi->ADH_COMM);
    error_vx=error_vx_max;
    error_vy=error_vy_max;
    error_h=error_h_max;
#endif
    if(sm->grid->smpi->myid==0){
        
        FILE *fp;
        fp = fopen("error.out", "w");
        fprintf(fp,"x-velocity abs error: %30.20e :: relative error: %30.20e \n",error_vx/sm->grid->macro_nnodes,error_vx/sm->grid->macro_nnodes/analytic_vx);
        fprintf(fp,"y-velocity abs error: %30.20e :: relative error: %30.20e \n",error_vy/sm->grid->macro_nnodes,error_vy/sm->grid->macro_nnodes/analytic_vy);
        fprintf(fp,"head abs error: %30.20e :: relative error: %30.20e \n",error_h/sm->grid->macro_nnodes,error_h/sm->grid->macro_nnodes/1.);
        if (error_vx/sm->grid->macro_nnodes > 1.e-6 || error_vx/sm->grid->macro_nnodes/analytic_vx > 1.e-6)
            fprintf(fp,"VelX FAIL\n");
        else
            fprintf(fp, "PASS\n");
        
        if (error_vy/sm->grid->macro_nnodes > 1.e-6 || error_vy/sm->grid->macro_nnodes/analytic_vy > 1.e-6)
            fprintf(fp,"VelY FAIL\n");
        else
            fprintf(fp, "PASS\n");
        
        fclose(fp);
    }
    
}

