/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  sw2_bc_discharge.c This file collections functions responsible for
 *          the 2D shallow water boundary contributions to the elemental residual.          */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
static int printFieldWidth = 20;
static int printPrecision  = 10;
static int DEBUG_LOCAL = OFF;
static int DEBUG_PICKETS = OFF;
static int DEBUG_NODE_ID = UNSET_INT;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes the 1D element residual for the SW-2D model.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  \note CJT \:: Uses the SW2D model struct
 *
 * @param[out] elem_rhs      the 1D elemental residual array
 * @param[in]  mod           a pointer to the model struct
 * @param[in]  ie            the elemental id
 * @param[in]  pertubation   the Newton pertubation
 * @param[in]  perturb_node  the node to be perturbed
 * @param[in]  perturb_var   the variable to be perturbed
 * @param[in]  perturb_sign  the direction of Newton perturbation
 * @param[in]  DEBUG         a debug flag
 *
 *  \details Integrates the weak, discrete outflow boundary terms: \n
 *  \f{eqnarray*}{
 *   \residDA{i}{}{c}   &=& dt * \intbe{\phidd{i} \, q_n}{e} = dt * \intbe{\phidd{i} \, (\velb{h} \cdot \nrml) \, \depth{h}}{e} \\
 *   \residDA{i}{}{mx}  &=& dt * \bigg[
 *                                \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \ub{h} \, \depth{h})}
 *                             +  \bcFriction{}{e}{\phidd{i}}{\ub{h}}{\velb{h}}
 *                             +  \bcPressureDD{}{e}{\phidd{i}}{(\depth{h})}{x}
 *                             +  \bcPressureDD{}{e}{\phidd{i} \,\rho^h}{(\depth{h})}{x}
 *                               \bigg] \\
 *   \residDA{i}{}{my}  &=& dt * \bigg[
 *                                \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \ub{h} \, \depth{h})}
 *                             +  \bcFriction{}{e}{\phidd{i}}{\vb{h}}{\velb{h}}
 *                             +  \bcPressureDD{}{e}{\phidd{i}}{(\depth{h})}{y}
 *                             +  \bcPressureDD{}{e}{\phidd{i} \,\rho^h}{(\depth{h})}{y}
 *                               \bigg]
 *  \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int fe_sw2_bc_discharge(SMODEL_SUPER *mod, double *elem_rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG) {
    int i;
#ifdef _DEBUG
    if (DEBUG_NODE_ID != UNSET_INT) {
        for (i=0; i<mod->grid->elem1d[ie].nnodes; i++) {
            if (mod->grid->elem1d[ie].nodes[i] == DEBUG_NODE_ID - 1)  {DEBUG = ON; break;} else {DEBUG = OFF; DEBUG_LOCAL = OFF; DEBUG_PICKETS = OFF;}
        }
    }
    if (DEBUG == ON) DEBUG_LOCAL = ON;
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__,__LINE__);
    time_t time1;  time(&time1);
#endif
    // aliases
    SSW *sw2 = mod->sw;
    SELEM_1D elem1d = mod->grid->elem1d[ie];
    // instead of string, this will come from physics mat or other coverage?
    //regardless this will always be in mat1d
    //int string = elem1d.string;
    int string = mod->elem1d_physics_mat[ie];
    int nnodes = elem1d.nnodes;
    int ndof = nnodes*3; //this is a 3 dof problem
    SVECT2D nrml = elem1d.nrml;
    double dt = mod->ts->dt;//need to pull this outside of loop
    double djac = elem1d.djac;
    double g = mod->gravity;
    STR_VALUE *str_values = mod->str_values;//what is this?
    //SWEIR_C *weir = mod->weir; //how should we handle these?
    //SFLAP_C *flap = mod->flap;
    //SSLUICE_C *sluice = mod->sluice; 
#ifdef _DEBUG
    assert(nnodes == NDONSEG); // only 1D element is a line segment
    assert(djac > 0); // only 1D element is a line segment
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // INDEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    if (perturb_var == PERTURB_NONE) perturb_sign = 0;
    double elem_head[nnodes];
    //global_to_local_dbl_cg(elem_head, mod->sol, elem1d.nodes, nnodes, PERTURB_H, mod->dof_map_local, mod->node_physics_mat);
    if (perturb_var == PERTURB_H) elem_head[perturb_node] += perturb_sign * perturbation;
    SVECT2D elem_vel[nnodes];
    //global_to_local_SVECT2D_cg(elem_vel, mod->sol, elem1d.nodes, nnodes, PERTURB_U, PERTURB_V, mod->dof_map_local, mod->node_physics_mat);
    if (perturb_var == PERTURB_U) {
        elem_vel[perturb_node].x += perturb_sign * perturbation;
    }
    if (perturb_var == PERTURB_V) {
        elem_vel[perturb_node].y += perturb_sign * perturbation;
    }
    double u_tan = one_2 * (elem_vel[0].x + elem_vel[1].x) * nrml.y - one_2 * (elem_vel[0].y + elem_vel[1].y) * nrml.x;
    double mag_utan = fabs(u_tan);
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEPENDENT VARIABLES
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    double u[NDONSEG], v[NDONSEG];
    for (i=0; i<nnodes; i++) { u[i] = elem_vel[i].x; v[i] = elem_vel[i].y; }
    double h_avg = one_2 * (elem_head[0] + elem_head[1]); if (h_avg < SMALL) h_avg = SMALL;
    double roughness[nnodes];
    roughness[0] = fe_sw2_get_roughness(mod, elem_head[0], 0.0, string, UNUSED);
    roughness[1] = fe_sw2_get_roughness(mod, elem_head[1], 0.0, string, UNUSED);
    double resistance[nnodes];
    resistance[0] = one_2 * roughness[0] * mag_utan;
    resistance[1] = one_2 * roughness[1] * mag_utan;
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // DEBUG SCREEN PRINT
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        printf("SW-2D 1D ELEM RESID :: ie: %d \t dt: %20.10f \t djac: %20.10f nrml: {%10.5e, %10.5e}\n",ie,dt,djac,nrml.x,nrml.y);
        if (perturb_var == PERTURB_H) {
            printf("\t perturbing H || node: %d || perturbation: %20.10e\n",elem1d.nodes[perturb_node],perturb_sign*perturbation);
        } else if (perturb_var == PERTURB_U) {
            printf("\t perturbing U || node: %d || perturbation: %20.10e\n",elem1d.nodes[perturb_node],perturb_sign*perturbation);
        } else if (perturb_var == PERTURB_V) {
            printf("\t perturbing V || node: %d || perturbation: %20.10e\n",elem1d.nodes[perturb_node],perturb_sign*perturbation);
        }
        printf("\n--------------------------------------------------------- \n");
        printScreen_debug2_dbl("elem_head", elem_head, nnodes, elem1d.nodes);
        printScreen_debug_svec2d("elem_vel", elem_vel, nnodes, elem1d.nodes);
    }
#endif
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                FINITE ELEMENT INTEGRATIONS
     *==========================================================================================*/
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                 DISCHARGE CONTRIBUTION
     *--------------------------------------------------------------------------------------------
     * Calculates the natural discharge addition to the elemental residual. \n
     * \note
     *********************************************************************************************/
    if (elem_head[0] < SMALL || elem_head[1] < SMALL)  return 0;
    double discharge[2] = { 0., 0. };
    double troughness[2] = { 0., 0. };
    double fact = 0., convey = 0., three_2 = 3.0 / 2.0;
    convey = str_values[string].conveyance;
    troughness[0] = fe_sw2_get_roughness(mod, elem_head[0], mag_utan, string, UNUSED);
    troughness[1] = fe_sw2_get_roughness(mod, elem_head[1], mag_utan, string, UNUSED);
    int isers = str_values[string].ol_flow.isigma;
    fact = -1.0 * sseries_get_value(isers, mod->series_head,0) * sqrt(2. * g / troughness[0]) / convey;
    discharge[0] = fact * pow(elem_head[0], three_2);
    fact = -1.0 * sseries_get_value(isers, mod->series_head,0) * sqrt(2. * g / troughness[1]) / convey;
    discharge[1] = fact * pow(elem_head[0], three_2);
    double elem_avg_discharge = one_2 * (discharge[0] + discharge[1]);
    fe_sw2_1d_explicit_flow(elem_rhs, elem1d, ie, nrml, djac, dt, elem_avg_discharge, h_avg, u, v, DEBUG, DEBUG_LOCAL);
    fe_sw2_1d_wall_friction(elem_rhs, elem1d, ie, nrml, djac, dt, u, v, resistance, DEBUG, DEBUG_LOCAL);
#ifdef _DEBUG
        if (DEBUG == ON || DEBUG_LOCAL == ON) printf("DEBUG BCT_DIS_NEU :: discharge: %20.10e, %20.10e \n",discharge[0],discharge[1]);
#endif
}
