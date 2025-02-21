/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  sdt.c This file collects methods of the SDT structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Intializes a SDT structure
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] t          (SDT *)  a pointer to an AdH design-level SDT structure
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void sdt_init(SDT *t) {
	t->t_init = 0.;
	t->t_prev = 0.;
	t->t_final = 0.;
	t->tau_temporal = 0.;
	t->dt = 0.;
	t->dt_old = 0.;
	t->dt_err = 0.;
	t->dt_prev = 0.;
	t->ientry_out = 0; 
	t->green_ampt = FALSE;
    t->time = 0.0;
    t->nt = 0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculated time-unit conversation factor
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] units (int) units to convert
 * @param[in] direct (int) TO or FROM seconds
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int sdt_get_conversion_factor(int units,int direct) {
	if (direct == TO) {
		switch (units) {
			case SECONDS:
				return 1.;
			case MINUTES:
				return 60.;
			case HOURS:
				return 3600.;
			case DAYS:
				return 86400.;
			case WEEKS:
				return 604800.;
			case MONTHS:
				return 2629744.;
			case YEARS:
				return 31556926.;
			default:
				return -1.;
		}
	} else if (direct == FROM) {
		switch (units) {
			case SECONDS:
				return 1;
			case MINUTES:
				return pow(60.0, -1.0);
			case HOURS:
				return pow(3600.0, -1.0);
			case DAYS:
				return pow(86400.0, -1.0);
			case WEEKS:
				return pow(604800.0, -1.0);
			case MONTHS:
				return pow(2629744.0, -1.0);
			case YEARS:
				return pow(31556926.0, -1.0);
			default:
				return -1;
		}
	} else {
		return -1;
	}
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Checks user input/output time units
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] uin  (int) input units
 * @param[in] uout (int) output units
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void sdt_check_units(int uin, int uout) {
  if (uin < SECONDS || uin > WEEKS) {
    //printf("\nSERIES %d ERROR\n", series.id);
    //sseries_printScreen(series, 0);
    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
    tl_error(">> Invalid input time units. Use a value between 0 and 4.");
  }
  if (uout < SECONDS || uout > WEEKS) {
    //printf("\nSERIES %d ERROR\n", series.id);
    //sseries_printScreen(series, 0);
    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
    tl_error(">> Invalid output time units. Use a value between 0 and 4.");
  }
}
