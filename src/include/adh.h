//#ifndef H_SADH_
//#define H_SADH_

void print_build_info();

/* standard header files */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <stdlib.h>
#include <stddef.h>
#include <limits.h>
#include <ctype.h>
#include <stdbool.h>
//Mark added
#include <umfpack.h>
#include <scotch.h>

#ifdef _HDF5
#include <hdf5.h>
#endif

#ifdef _MPI
#include <mpi.h>
#endif
#include "messg.h"

#include "vars.h"
#include "define.h"
#include "macro.h"
#include "model_codes.h"

#ifdef _PETSC
//Mark changed
#include <petsc.h>
//#include <petscksp.h>
//#include <petscts.h>
#endif

#include "debug.h"
#include "header_tl_alloc.h"

#include "assert.h"
//#include "constants.h"

// STRUCTURES
#include "sdt.h"
#include "scoverage.h"
#include "svect2d.h"
#include "svect.h"
#include "snode.h"
#include "stensor.h"
#include "selem_1d.h"
#include "selem_2d.h"
#include "selem_3d.h"
#include "squad.h"
#include "smodel.h"
#include "slist_items.h"
#include "smpi.h"
#include "smeteor.h"
#include "sflags.h"
#include "sstr_value.h"
#include "sgrid.h"
#include "sarray.h"
#include "tokens.h"
#include "sivar_position.h"
#include "smat_grid.h"
#include "smat_sw.h"
#include "smat_gw.h"
#include "smat_transport.h"
#include "smat_physics.h"
#include "slin_sys.h"
#include "dofmaps.h"
#include "sdvar.h"
#include "ssw.h"
#include "sfile.h"
#include "sseries.h"
#include "smodel.h"
#include "smodel_super.h"
#include "smodel_design.h"


//Mark added physics module
#include "fe.h"
#include "sw2.h"
#include "poisson.h"
#include "heat.h"
#include "no_model.h"
#include "models.h"

#include "bc.h"

//Mark added
#include "residual.h"
#include "jacobian.h"
#include "la.h"
#include "newton.h"

//Mark added time loop
#include "time_loop.h"


// FOLDERS

//Mark, does test need to be moved up?
#include "tools.h"
#include "fnctn_xdmf.h"
#include "fr_defs.h"

//Mark added
#include "test.h"

//Adding function pointers that will always be needed
int (*adh_resid_routines[N_RESID_ROUTINES])(SMODEL_SUPER *, double *, int, double, int, int, int, int);
int (*adh_init_routines[N_INIT_ROUTINES])(SMODEL_SUPER *);
int (*adh_time_stepper[N_TIME_STEPPERS]) (SMODEL_SUPER*);
//#endif
