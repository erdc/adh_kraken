/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  scoverage.c This file collects methods of the SCOVERAGE structure for coverages   */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = ON;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates and intializes a SCOVERAGE structure
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] c (SCOVERAGE **)  a double pointer to an AdH SCOVERAGE structure
 * @param[in]  nelems1D (int) total # of 1D elements on the designer grid
 * @param[in]  nelems2D (int) total # of 2D elements on the designer grid
 * @param[in]  nelems3D (int) total # of 3D elements on the designer grid
 * 
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void scoverage_alloc_init(SCOVERAGE **coverage,int nelems1D,int nelems2D,int nelems3D, 
	int ncoverages1D,int ncoverages2D,int ncoverages3D) {

	// check input
	if (ncoverages1D + ncoverages2D + ncoverages3D < 1) return;
	if (ncoverages1D > 0 && nelems1D < 1) {tl_error("There are 1D coverages but no 1D elements.\n");}
	if (ncoverages2D > 0 && nelems2D < 1) {tl_error("There are 2D coverages but no 2D elements.\n");}
	if (ncoverages3D > 0 && nelems3D < 1) {tl_error("There are 3D coverages but no 3D elements.\n");}

	(*coverage) = (SCOVERAGE *) tl_alloc(sizeof(SCOVERAGE), 1);
	SCOVERAGE *c = (*coverage);

	c->nelems1D = nelems1D;
	c->nelems2D = nelems2D;
	c->nelems3D = nelems3D;
	c->ncoverages1D = ncoverages1D;
	c->ncoverages2D = ncoverages2D;
	c->ncoverages3D = ncoverages3D;

	if (c->ncoverages1D > 0) {c->coverage_1D = allocate_dptr_int(c->ncoverages1D,c->nelems1D);}
	if (c->ncoverages2D > 0) {c->coverage_2D = allocate_dptr_int(c->ncoverages2D,c->nelems2D);}
	if (c->ncoverages3D > 0) {c->coverage_3D = allocate_dptr_int(c->ncoverages3D,c->nelems3D);}
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees a SCOVERAGE structure
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] c (SCOVERAGE *)  a pointer to an AdH SCOVERAGE structure
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void scoverage_free(SCOVERAGE *c) {
	free_dptr_int(c->coverage_1D,c->nelems1D,c->ncoverages1D);
	free_dptr_int(c->coverage_2D,c->nelems2D,c->ncoverages2D);
	free_dptr_int(c->coverage_3D,c->nelems3D,c->ncoverages3D);
	c = (SCOVERAGE *) tl_free(sizeof(SCOVERAGE), 1, c);
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads coverages from file
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] c (SCOVERAGE *)  a pointer to an AdH SCOVERAGE structure
 * @param[in]  fp (FILE *) a pointer to a coverage file
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void scoverage_read(SCOVERAGE **cover, char *filebase) {
	int i;
	char *line = NULL;
	char *token;
	size_t len = 0;
	ssize_t read;

	SFILE file;
    const char suffix[6] = ".cover", mode[1] = "r";
    sfile_open(&file,filebase,NULL,NULL,suffix,mode,TRUE);
    FILE *fp = file.fp;

	// read for counts (just body elements?)
	int n1d=0, n2d=0, n3d=0, ne1d=0, ne2d=0, ne3d=0;
	while ((read = getline(&line,&len,fp)) != -1) {
		get_token(line,&token); if (token == NULL) continue;
		if (strcmp(token,"NCOVERAGES") == 0) {
			n1d = get_next_token_int(&token);
			n2d = get_next_token_int(&token);
			n3d = get_next_token_int(&token);
		}
		if (strcmp(token,"ELEMS_SEG") == 0) {
			ne1d += get_next_token_int(&token);
		}
		if (strcmp(token,"ELEMS_TRI") == 0) {
			ne2d += get_next_token_int(&token);
		}
		if (strcmp(token,"ELEMS_QUAD") == 0) {
			ne2d += get_next_token_int(&token);
		}
		if (strcmp(token,"ELEMS_TET") == 0) {
			ne3d += get_next_token_int(&token);
		}
		if (strcmp(token,"ELEMS_PRISM") == 0) {
			ne3d += get_next_token_int(&token);
		}		
	}
	rewind(fp);


	// allocate
	scoverage_alloc_init(cover,ne1d,ne2d,ne3d,n1d,n2d,n3d);
	SCOVERAGE *c = *cover;

	// read to store
	n1d=0;n2d=0;n3d=0;
	while ((read = getline(&line, &len, fp)) != -1) {
		get_token(line,&token); if (token == NULL) continue;
		if (strcmp(token,"SEG") == 0) {
			for (i=0; i<c->ncoverages1D; i++) {
				c->coverage_1D[i][n1d] = get_next_token_int(&token);
			}
			n1d++;
		}
		if (strcmp(token,"TRI") == 0) {
			for (i=0; i<c->ncoverages2D; i++) {
				c->coverage_2D[i][n2d] = get_next_token_int(&token);
			}
			n2d++;
		}
		if (strcmp(token,"QUAD") == 0) {
			for (i=0; i<c->ncoverages2D; i++) {
				c->coverage_2D[i][n2d] = get_next_token_int(&token);
			}
			n2d++;
		}
		if (strcmp(token,"TET") == 0) {
			for (i=0; i<c->ncoverages3D; i++) {
				c->coverage_3D[i][n3d] = get_next_token_int(&token);
			}
			n3d++;
		}
		if (strcmp(token,"PRISM") == 0) {
			for (i=0; i<c->ncoverages3D; i++) {
				c->coverage_3D[i][n3d] = get_next_token_int(&token);
			}
			n3d++;
		}	
	}
	if (c->nelems1D != n1d) {tl_error("Error with 1D element listing in coverage file!\n");}
	if (c->nelems2D != n2d) {tl_error("Error with 2D element listing in coverage file!\n");}
	if (c->nelems3D != n3d) {tl_error("Error with 3D element listing in coverage file!\n");}

	if (DEBUG){
		printf("Params File Info:\n");
		printf("-- 1D element count: %10d || 1D coverages: %d \n",c->nelems1D,c->ncoverages1D);
		printf("-- 2D element count: %10d || 2D coverages: %d \n",c->nelems2D,c->ncoverages2D);
		printf("-- 3D element count: %10d || 3D coverages: %d \n",c->nelems3D,c->ncoverages3D);
	}

	fclose(fp);

}