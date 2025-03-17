/* This routine initializes dt for the simulation */
#include "adh.h"

void tc_init(SMODEL_DESIGN *mod)
{
  int interval;                 /* the interval of the new time step */

  mod->ts.dt_old = mod->ts.dt;

  //need to fix this, Mark

  /* initializes the new dt to the time series provided by the user */
  interval = sseries_get_interval(*(mod->series_dt), mod->ts.t_prev);
  mod->ts.dt = tc_eval_series(*(mod->series_dt), interval, mod->ts.t_prev, 0);//

  /* if we are within one time step of the end, then take exactly the right amount */
  if (mod->ts.t_prev + mod->ts.dt >= mod->ts.t_final)
    mod->ts.dt = mod->ts.t_final - mod->ts.t_prev;//

  /* if we are within one and a little time step of the end, then take half the time */
  else if (mod->ts.t_prev + 2.0 * mod->ts.dt >= mod->ts.t_final)
    mod->ts.dt = 0.5 * (mod->ts.t_final - mod->ts.t_prev);//

  /* check for a ridiculous time step */
  if (mod->ts.dt < SMALL) {
      printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
      tl_error("Time step is too small.");
  }

  /* makes sure everyone is on the same dt */
  //old_dt = messg_dmin(old_dt);
  //dt = messg_dmin(dt);
}
