/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  This fil initializes and SMODEL_DESIGN structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = ON;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads an AdH Designer Model
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] dmod     (SMODEL_DESIGN **) a pointer to an AdH design-level model
 * @param[in]    filename (CHAR *) the design model filename

 * \note CJT :: Assume all model parameters within a SuperModel are given in the same file
 * \ note: CJT:: New AdH Input Files - .bc, .init, .cov, .geo
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void smodel_design_read(SMODEL_DESIGN *dmod, char *filename) {
    
    char *ext = NULL; 
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *token;
    char type[40], fn[40];
    int imono=0,model1=0,model2=0,nmats=0,imat=UNSET_INT,ntrns=0,ie=0;
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Open Designer File
    //++++++++++++++++++++++++++++++++++++++++++++++
    SFILE file;
    char mode[1] = "r";
    sfile_open(&file,filename,NULL,NULL,NULL,mode,TRUE);
    FILE *fp = file.fp;

    //++++++++++++++++++++++++++++++++++++++++++++++
    // Read DEBUG/SCREEN/FILE OUTPUT options
    //++++++++++++++++++++++++++++++++++++++++++++++
    while ((read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token);
        if (strcmp(token, "DEBUG") == 0) {
            //read_bc_DEBUG(sm,&token);
        } else if (strcmp(token, "SCREEN_OUTPUT") == 0) {
            //read_bc_SCREEN_OUTPUT(sm,&token);
        }
    }
    rewind(fp);
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Read Grid Filename
    //++++++++++++++++++++++++++++++++++++++++++++++
    while ((read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token); //printf("token: %s\n",token);
        if (token == NULL) continue;
        if (strcmp(token,"GRID") == 0) {
            get_next_token(&token);
            strcpy(dmod->filename_grid, token);
            if (DEBUG) printf("filename_grid found: %s\n",dmod->filename_grid);
        }
    }
    rewind(fp);

    //++++++++++++++++++++++++++++++++++++++++++++++
    // Read Parameter Coverage Filename (SuperModel Shared)
    //++++++++++++++++++++++++++++++++++++++++++++++
    while ((read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token); //printf("token: %s\n",token);
        if (token == NULL) continue;
        if (strcmp(token,"COVERAGE") == 0) {
            get_next_token(&token);
            strcpy(dmod->filename_params, token);
        }
    }
    rewind(fp);

    //++++++++++++++++++++++++++++++++++++++++++++++
    // Read Time and Output Information
    //++++++++++++++++++++++++++++++++++++++++++++++
    int in_units, out_units, ientry;
    while ((read = getline(&line, &len, fp)) != -1) {
        if (get_compare_token(line,&token,"TC") == 1) {
            get_next_token(&token); if (token == NULL) continue;
            if (strcmp(token,"T0") == 0) {
                dmod->ts.t_init = get_next_token_int(&token); //printf("token: %s\n",token);
                in_units = get_next_token_int(&token); //printf("token: %s\n",token);
                dmod->ts.t_init *= sdt_get_conversion_factor(in_units,TO);
                dmod->ts.t_prev = dmod->ts.t_init;
            } else if (strcmp(token,"TF") == 0) {
                dmod->ts.t_final = get_next_token_int(&token); //printf("token: %s\n",token);
                in_units = get_next_token_int(&token); //printf("token: %s\n",token);
                dmod->ts.t_final *= sdt_get_conversion_factor(in_units,TO);
            }
        }
        if (get_compare_token(line,&token,"SERIES") == 1) {
            get_next_token(&token); if (token == NULL) continue;
            printf("token: %s\n",token);
            if (strcmp(token,"DT") == 0) {
                dmod->series_dt = sseries_read_allocate(fp,NULL,NULL,DT_SERIES,&token,UNSET_INT); 
                //dmod->series_dt = sseries_read_allocate(NULL,NULL,DT_SERIES,&token,UNSET_INT); 
            } else if (strcmp(token,"WRITE") == 0) {
                
#ifdef _MESSG
                tag(dmod->grid->smpi->ADH_COMM);
#else
                tag();
#endif
                //dmod->series_out = sseries_read_allocate(NULL,NULL,OUTPUT_SERIES,&token,UNSET_INT); 
                dmod->series_out = sseries_read_allocate(fp,NULL,NULL,OUTPUT_SERIES,&token,UNSET_INT); 
#ifdef _MESSG
                tag(dmod->grid->smpi->ADH_COMM);
#else
                tag();
#endif
            } else if (strcmp(token,"AWRITE") == 0) {
                //dmod->series_out = sseries_read_allocate(NULL, NULL,OUTPUT_SERIES,&token,UNSET_INT); 
            }
        }
    }
    rewind(fp);

    //++++++++++++++++++++++++++++++++++++++++++++++
    // Count SuperModels
    //++++++++++++++++++++++++++++++++++++++++++++++
    dmod->nSuperModels = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        if (get_compare_token(line,&token,"MONO") == 1) {
            imono = get_next_token_int(&token);
            if (imono > dmod->nSuperModels) dmod->nSuperModels = imono;
        }
    }
    printf("%d Super models found\n",dmod->nSuperModels);
    rewind(fp);
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Allocate and Initialize SuperModels
    //++++++++++++++++++++++++++++++++++++++++++++++
    smodel_super_alloc_init(&dmod->superModel,dmod->nSuperModels);
    if (dmod->nSuperModels > 1) {
        dmod->superCouple = (SMODEL_COUPLE *) tl_alloc(sizeof(SMODEL_COUPLE), dmod->nSuperModels-1);
    }
        
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Read SuperModel base filenames
    //++++++++++++++++++++++++++++++++++++++++++++++
    while ((read = getline(&line, &len, fp)) != -1) {
        if (get_compare_token(line,&token,"MONO") == 1) {
            imono = get_next_token_int(&token) - 1; //printf("token: %s\n",token);
            dmod->superModel[imono].id = imono;
            get_next_token(&token); //printf("token: %s\n",token);
            strcpy(dmod->superModel[imono].filebase, token);
            printf("SuperModel: %d || Filename: %s\n",imono,dmod->superModel[imono].filebase);
        }
    }
    rewind(fp);
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Read SuperModel Coupling Types
    //++++++++++++++++++++++++++++++++++++++++++++++
    int icpl = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token); //printf("token: %s\n",token);
        if (token == NULL) continue;
        if (strcmp(token,"COUPLE") == 0) {
            dmod->superCouple[icpl].model1 = get_next_token_int(&token);
            dmod->superCouple[icpl].model2 = get_next_token_int(&token);
            get_next_token(&token);
            if (strcmp(token,"LAG") == 0) {
                dmod->superCouple[icpl].type = COUPLING_LAG;
            } else if (strcmp(token,"FLUX") == 0) {
                dmod->superCouple[icpl].type = COUPLING_FLUX;
            } else  if (strcmp(token,"STRANG") == 0) {
                dmod->superCouple[icpl].type = COUPLING_STRANG;
            }
            printf("Coupling Model %d to %d Using: %d\n",
                   dmod->superCouple[icpl].model1,
                   dmod->superCouple[icpl].model2,
                   dmod->superCouple[icpl].type);
            icpl++;
        }
    }
    rewind(fp);

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    
    fclose(fp); fp = NULL;

    //++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++
    return;
}
