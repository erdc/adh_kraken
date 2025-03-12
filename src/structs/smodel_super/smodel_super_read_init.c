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
 * \note: CJT :: Must be called after BC read to setup which variables should be read
 * \note: CJT :: DataSet has a structured format
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void write_stats(int ndim, double *min, double *max);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void smodel_super_read_init(SMODEL_SUPER *sm, char *filebase) {
    
    int i,ivar,ndim,ndiml,ip,inode,lnnodes,lnelems;
    char *token,mode[MAXLINE],ext[MAXLINE],name[MAXLINE],preName[MAXLINE],str[MAXLINE];
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    double max[3], min[3];
    bool found = false;
    int nnodes = sm->grid->nnodes;
    SIVAR_POSITION *ivar_pos = &sm->ivar_pos;
    int **ivars = sm->ivars;
    double *s = NULL;
    SDVAR_POSITION *dvar_pos = &sm->dvars.sdvar_pos_node;
    double **dvars = sm->dvars.nodal_dvar;
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Open Initialization File
    //++++++++++++++++++++++++++++++++++++++++++++++
    SFILE file;
    strcpy(mode,"r"); strcpy(ext,".init");
    sfile_open(&file,filebase,NULL,NULL,ext,mode,TRUE);
    FILE *fp = file.fp; // alias
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Read INIT file
    //++++++++++++++++++++++++++++++++++++++++++++++
    while ((read = getline(&line, &len, fp)) != -1) { //printf("line: %s\n",line);
        get_token(line,&token); if (token == NULL) continue; //printf("token: %s\n",token);
        if (strcmp(token,"DATASET") == 0) {
            found = false;
            strcpy(name,"");
            
            // ++++++++++++++++++++++++++
            // Read DataSet Header Info
            // ++++++++++++++++++++++++++
            read = getline(&line, &len, fp); get_token(line,&token); assert(strcmp(token,"OBJTYPE") == 0);
            read = getline(&line, &len, fp); get_token(line,&token); assert(strcmp(token,"ND") == 0);
            lnnodes = get_next_token_int(&token); assert(lnnodes = nnodes);
            read = getline(&line, &len, fp); get_token(line,&token); assert(strcmp(token,"NC") == 0);
            read = getline(&line, &len, fp); get_token(line,&token); assert(strcmp(token,"DIM") == 0);
            ndim = get_next_token_int(&token); assert(ndim > 0 && ndim < 3);
            read = getline(&line, &len, fp); get_token(line,&token); assert(strcmp(token,"NAME") == 0);
            get_next_token(&token);
            if (strcmp(token,"OLD") == 0) {
                strcpy(preName,token);
                s = sm->sol_old;
                get_next_token(&token); if (token == NULL) continue;
            } else if (strcmp(token,"OLDER") == 0) {
                strcpy(preName,token);
                s = sm->sol_older;
                get_next_token(&token); if (token == NULL) continue;
            } else {
                strcpy(preName,"");
                s = sm->sol;
            }
            strcpy(name,token);
            read = getline(&line, &len, fp); get_token(line,&token); assert(strcmp(token,"TS") == 0);
            
            //if (DEBUG == ON) {printf("<<< reading DATASET for variable: %s %s of dim: %d\n",preName,name,ndim);}
            
            // +++++++++++++++++++++++++++
            // Independent Variables Read
            // Stores into solution array
            // +++++++++++++++++++++++++++
            for (ivar=0; ivar<N_IVARS_TOTAL; ivar++) { //printf("var_name: %s\n",ivar_pos->var_name[ivar]);
                //printf("name: %s || IVAR_NAME[ivar]: %s\n",name,IVAR_NAME[ivar]);
                if (strcmp(name,IVAR_NAME[ivar]) == 0) {
                    if (DEBUG) printf("---- initializing: %s %s\n",preName,IVAR_NAME[ivar]);
                    
                    for (i=0; i<3; i++) {max[i] = -99999999.; min[i] = 99999999.;}
                    for (inode=0; inode<nnodes; inode++) {
                        read = getline(&line, &len, fp); //printf("line: %s\n",line); //exit(-1);
                        get_token(line,&token);
                        ndiml = 0;
                        while (token != NULL) {
                            ip = ivar_pos->var[ivar + ndiml];
                            if (ivars[ip][inode] == UNSET_INT) continue;
                            sscanf(token, "%lf", &s[ivars[ip][inode]]); // independent variables go into the solution array
                            if (s[ivars[ip][inode]] > max[ndiml]) {max[ndiml] = s[ivars[ip][inode]];}
                            if (s[ivars[ip][inode]] < min[ndiml]) {min[ndiml] = s[ivars[ip][inode]];}
                            ndiml++;
                            get_next_token(&token);
                        }
                        assert(ndiml == ndim);
                    }
                    read = getline(&line, &len, fp); get_token(line,&token); assert(strcmp(token,"ENDDS") == 0);
                    if (DEBUG == ON) {write_stats(ndim,min,max);}
                    found = true; break;
                }
            }
            if (found) {
                //if (DEBUG == ON) {printf("<<< Finished reading DATASET || variable: %s %s || dim: %d \n",preName,name,ndim);}
                continue;
            }
            
            // +++++++++++++++++++++++++++
            // Dependent Variables Read
            // Stores into dvar matrix
            // +++++++++++++++++++++++++++
            for (ivar=0; ivar<N_DVARS; ivar++) { //printf("var_name: %s\n",ivar_pos->var_name[ivar]);
                if (strcmp(name,DVAR_NAME[ivar]) == 0) {
                    if (DEBUG) printf("---- initializing: %s %s\n",preName,DVAR_NAME[ivar]);
                    
                    for (i=0; i<3; i++) {max[i] = -99999999.; min[i] = 99999999.;}
                    for (inode=0; inode<nnodes; inode++) {
                        read = getline(&line, &len, fp); //printf("line: %s\n",line); //exit(-1);
                        get_token(line,&token);
                        ndiml = 0;
                        while (token != NULL) {
                            ip = dvar_pos->var[ivar + ndiml];
                            //printf("token: %s || ivar: %d || ip: %d || ndiml: %d\n",token,ivar + ndiml,ip,ndiml);
                            if (ip == UNSET_INT) continue;
                            sscanf(token, "%lf", &dvars[ip][inode]); // dependent variables are stored in a matrix
                            if (dvars[ip][inode] > max[ndiml]) {max[ndiml] = dvars[ip][inode];}
                            if (dvars[ip][inode] < min[ndiml]) {min[ndiml] = dvars[ip][inode];}
                            ndiml++;
                            get_next_token(&token);
                        }
                        assert(ndiml == ndim);
                    }
                    read = getline(&line, &len, fp); get_token(line,&token); assert(strcmp(token,"ENDDS") == 0);
                    write_stats(ndim,min,max);
                    found = true; break;
                }
            }
            if (!found) {
                sprintf(str, "WARNING: Initialization variable %s %s not found by AdH.\n",preName,name);
                tl_error(str);
            } else {
                if (DEBUG == ON) {printf("<<< Finished reading DATASET || variable: %s %s || dim: %d \n",preName,name,ndim);}
            }
        }
    }
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    //tl_check_all_pickets(__FILE__,__LINE__);
    //exit(1);
    
    fclose(file.fp);
    
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void write_stats(int ndim, double *min, double *max) {
    for (int i=0; i<ndim; i++) {
        printf("---- dim[%d] || min/max: %f/%f\n",i,min[i],max[i]);
    }
}
