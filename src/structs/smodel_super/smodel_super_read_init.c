/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  For reading superModel initializations */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = ON;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads an AdH SuperModel initialization file
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] sm     (SMODEL_DESIGN **) a pointer to an AdH design-level model
 * @param[in]    filename (CHAR *) the design model filename
 * 
 * \note: CJT :: New AdH Input Files - .bc, .init, .cov, .geo
 * \note: CJT :: Need to figure out how to read dependent variables!
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_1col(double *sol, int **ivars, int nnodes, char **line, int varcode, FILE *fp, size_t *len) {
	ssize_t read;
	for (int inode=0; inode<nnodes; inode++) {
		read = getline(line, len, fp);
        //printf("varcode: %d || inode: %d ||  ivars[varcode][inode]: %d\n",varcode,inode,ivars[varcode][inode]);   
        if (ivars[varcode][inode] == UNSET_INT) continue;
        //printf("varcode: %d || inode: %d ||  ivars[varcode][inode]: %d\n",varcode,inode,ivars[varcode][inode]);   
        sscanf(line, "%lf", &sol[ivars[varcode][inode]]);
	}
}
void read_2col(double *sol, int **ivars, int nnodes, char **line, int varcode1, int varcode2, FILE *fp, size_t *len) {
	ssize_t read;
	for (int inode=0; inode<nnodes; inode++) {
		read = getline(line, len, fp);
        if (ivars[varcode1][inode] == UNSET_INT || 
            ivars[varcode2][inode] == UNSET_INT) continue;
        //printf("inode: %d || varcodes: [%d,%d] || ivars1[varcode1][inode]: %d || ivars2[varcode2][inode]: %d\n",
        //            inode,varcode1,varcode2,ivars[varcode1][inode],ivars[varcode2][inode]); 

        sscanf(line, "%lf %lf", &sol[ivars[varcode1][inode]],&sol[ivars[varcode2][inode]]);
	}
}
void read_3col(double *sol, int **ivars, int nnodes, char **line, int varcode1, int varcode2, int varcode3, FILE *fp, size_t *len) {
	ssize_t read;
	for (int inode=0; inode<nnodes; inode++) {
        if (ivars[varcode1][inode] == UNSET_INT || 
            ivars[varcode2][inode] == UNSET_INT ||
            ivars[varcode3][inode] == UNSET_INT) continue;
		//read = getline(line, len, fp);
		sscanf(line, "%lf %lf %lf", &sol[ivars[varcode1][inode]],&sol[ivars[varcode2][inode]],&sol[ivars[varcode3][inode]]);
	}
}

void smodel_super_read_init(SMODEL_SUPER *sm, char *filebase) {

    char *token,mode[MAXLINE],ext[MAXLINE];
    SFILE file;
    strcpy(mode,"r"); strcpy(ext,".init");
    sfile_open(&file,filebase,NULL,NULL,ext,mode,TRUE);
    FILE *fp = file.fp;

    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    double dtmp;
    int conID = UNSET_INT;

    int nnodes = sm->grid->nnodes; // alias
    SIVAR_POSITION *ivar_pos = &sm->ivar_pos; 
    double *s = sm->sol;
    double *s_old = sm->sol_old;
    double *s_older = sm->sol_older;   

    //++++++++++++++++++++++++++++++++++++++++++++++
    // Read INIT file
    //++++++++++++++++++++++++++++++++++++++++++++++
    while ((read = getline(&line, &len, fp)) != -1) {
        //printf("line: %s\n",line);
        get_token(line,&token); //printf("token: %s\n",token);
        if (token == NULL) continue;
        // //printf("token: %s\n",token);
        if (strcmp(token,"DEPTH") == 0) {
            if (DEBUG) printf("-- initializing depth\n");
        	read_1col(s,sm->ivars,nnodes,&line,ivar_pos->h,fp,&len);
        }
        if (strcmp(token,"OLD_DEPTH") == 0) {
            if (DEBUG) printf("-- initializing old depth\n");
        	read_1col(s_old,sm->ivars,nnodes,&line,ivar_pos->h,fp,&len);
        }
        if (strcmp(token,"OLDER_DEPTH") == 0) {
            if (DEBUG) printf("-- initializing older depth\n");
        	read_1col(s_older,sm->ivars,nnodes,&line,ivar_pos->h,fp,&len);
        }
        if (strcmp(token,"VELOCITY") == 0) {
            if (DEBUG) printf("-- initializing velocity\n");
        	read_3col(s,sm->ivars,nnodes,&line,ivar_pos->u,ivar_pos->v,ivar_pos->w,fp,&len);
        }
        if (strcmp(token,"OLD_VELOCITY") == 0) {
            if (DEBUG) printf("-- initializing old velocity\n");
        	read_3col(s_old,sm->ivars,nnodes,&line,ivar_pos->u,ivar_pos->v,ivar_pos->w,fp,&len);
        }
        if (strcmp(token,"OLDER_VELOCITY") == 0) {
            if (DEBUG) printf("-- initializing older velocity\n");
        	read_3col(s_older,sm->ivars,nnodes,&line,ivar_pos->u,ivar_pos->v,ivar_pos->w,fp,&len);
        }
        if (strcmp(token,"DA_VELOCITY") == 0) {
            if (DEBUG) printf("-- initializing depth-averaged velocity\n");
        	read_2col(s,sm->ivars,nnodes,&line,ivar_pos->uda,ivar_pos->vda,fp,&len);
        }
        if (strcmp(token,"OLD_DA_VELOCITY") == 0) {
            if (DEBUG) printf("-- initializing old depth-averaged velocity\n");
        	read_2col(s_old,sm->ivars,nnodes,&line,ivar_pos->uda,ivar_pos->vda,fp,&len);
        }
        if (strcmp(token,"OLDER_DA_VELOCITY") == 0) {
            if (DEBUG) printf("-- initializing older depth-averaged velocity\n");
        	read_2col(s_older,sm->ivars,nnodes,&line,ivar_pos->uda,ivar_pos->vda,fp,&len);
        }
        if (strcmp(token,"DISPLACEMENT") == 0) {
            if (DEBUG) printf("-- initializing displacement\n");
        	read_1col(s,sm->ivars,nnodes,&line,ivar_pos->dpl,fp,&len);
        }
        if (strcmp(token,"OLD_DISPLACEMENT") == 0) {
            if (DEBUG) printf("-- initializing old displacement\n");
        	read_1col(s_old,sm->ivars,nnodes,&line,ivar_pos->dpl,fp,&len);
        }
        if (strcmp(token,"OLDER_DISPLACEMENT") == 0) {
            if (DEBUG) printf("-- initializing older displacement\n");
        	read_1col(s_older,sm->ivars,nnodes,&line,ivar_pos->dpl,fp,&len);
        }
        if (strcmp(token,"PRESSURE") == 0) {
            if (DEBUG) printf("-- initializing pressure\n");
        	read_1col(s,sm->ivars,nnodes,&line,ivar_pos->prs,fp,&len);
        }
        if (strcmp(token,"OLD_PRESSURE") == 0) {
            if (DEBUG) printf("-- initializing old pressure\n");
        	read_1col(s_old,sm->ivars,nnodes,&line,ivar_pos->prs,fp,&len);
        }
        if (strcmp(token,"OLDER_PRESSURE") == 0) {
            if (DEBUG) printf("-- initializing older pressure\n");
        	read_1col(s_older,sm->ivars,nnodes,&line,ivar_pos->prs,fp,&len);
        }
        if (strcmp(token,"HEAT") == 0) {
            if (DEBUG) printf("-- initializing heat\n");
        	read_1col(s,sm->ivars,nnodes,&line,ivar_pos->dpl,fp,&len);
        }
        if (strcmp(token,"OLD_HEAT") == 0) {
            if (DEBUG) printf("-- initializing old heat\n");
        	read_1col(s_old,sm->ivars,nnodes,&line,ivar_pos->dpl,fp,&len);
        }
        if (strcmp(token,"OLDER_HEAT") == 0) {
            if (DEBUG) printf("-- initializing older heat\n");
        	read_1col(s_older,sm->ivars,nnodes,&line,ivar_pos->dpl,fp,&len);
        }
        if (strcmp(token,"SALINITY") == 0) {
            if (DEBUG) printf("-- initializing salinity\n");
        	read_1col(s,sm->ivars,nnodes,&line,ivar_pos->sal,fp,&len);
        }
        if (strcmp(token,"OLD_SALINITY") == 0) {
            if (DEBUG) printf("-- initializing old salinity\n");
        	read_1col(s_old,sm->ivars,nnodes,&line,ivar_pos->sal,fp,&len);
        }
        if (strcmp(token,"OLDER_SALINITY") == 0) {
            if (DEBUG) printf("-- initializing older salinity\n");
        	read_1col(s_older,sm->ivars,nnodes,&line,ivar_pos->sal,fp,&len);
        }
        if (strcmp(token,"CON") == 0) {
        	conID = get_next_token_int(&token);
            if (DEBUG) printf("-- initializing constituent #%d\n",conID);
        	read_1col(s,sm->ivars,nnodes,&line,ivar_pos->con[conID-1],fp,&len);
        }
        if (strcmp(token,"OLD_CON") == 0) {
        	conID = get_next_token_int(&token);
            if (DEBUG) printf("-- initializing old constituent #%d\n",conID);
        	read_1col(s_old,sm->ivars,nnodes,&line,ivar_pos->con[conID-1],fp,&len);
        }
        if (strcmp(token,"OLDER_CON") == 0) {
        	conID = get_next_token_int(&token);
            if (DEBUG) printf("-- initializing older constituent #%d\n",conID);
        	read_1col(s_older,sm->ivars,nnodes,&line,ivar_pos->con[conID-1],fp,&len);
        } 
    }

    //tl_check_all_pickets(__FILE__,__LINE__);
    //exit(1);

    fclose(file.fp);

}
