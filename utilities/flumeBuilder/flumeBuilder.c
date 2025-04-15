//#include "extrusion.h"  // includes AdH global header

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

char* concat(const char *s1, const char *s2) {
    char *result = malloc(strlen(s1) + strlen(s2) + 1); // +1 for the null-terminator
    // in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

typedef struct {
      double x, y, z;       /* the coordinates of the vector */
} SVECT;

typedef struct {
    int nMods;
    char model[10][50];
    double model_xmin[10];
    double model_xmax[10];
    double model_ymin[10];
    double model_ymax[10];
} SUPERMODEL;

SVECT svect_cross(SVECT v1, SVECT v2) {
    SVECT result;
    result.x = v1.y*v2.z - v1.z*v2.y;
    result.y = v1.z*v2.x - v1.x*v2.z;
    result.z = v1.x*v2.y - v1.y*v2.x;
    return result;
}

double svect_dotp(SVECT v1, SVECT v2) {
    return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}


int main(int argc, char *argv[]) {
    
    int i,j,k,igrp,iSupMod;

    if (strcmp(argv[1], "-help") == 0) {
        printf("usage: [root_file_name] [xmin] [xmax] [ymin] [ymax] [npx] [npy] [theta] [dz] [a] [b] [c] [flag_3D]\n");
        exit(0);
    }
    
    SVECT ref_vector;
    ref_vector.x = 0 - 0;
    ref_vector.y = 0 - 0;
    ref_vector.z = -1 - 0;
    
    char *file_name = argv[1];
    double xmin = atof(argv[2]);
    double xmax = atof(argv[3]);
    double ymin = atof(argv[4]);
    double ymax = atof(argv[5]);
    
    
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // CJT -- added for multi-physics building (must be hardwired for now)
    // DOESN'T WORK FOR ROTATIONS YET!
    int nSupMods = 2;
    SUPERMODEL sm[nSupMods];
    
    // SuperModel #1 Details
    sm[0].nMods = 2;
    strcpy(sm[0].model[0],"SW2");
    strcpy(sm[0].model[1],"DW");
    sm[0].model_xmin[0] = xmin;
    sm[0].model_xmax[0] = xmax/2.0;
    sm[0].model_ymin[0] = ymin;
    sm[0].model_ymax[0] = ymax;
    sm[0].model_xmin[1] = sm[0].model_xmax[0] + 1e-10;
    sm[0].model_xmax[1] = xmax;
    sm[0].model_ymin[1] = ymin;
    sm[0].model_ymax[1] = ymax;
    double h0 = 10;
    double u0 = 0;
    double v0 = 0;
    
    // SuperModel #2 Details
    sm[1].nMods = 1;
    strcpy(sm[1].model[0],"TRNS");
    sm[1].model_xmin[0] = xmin;
    sm[1].model_xmax[0] = xmax;
    sm[1].model_ymin[0] = ymin;
    sm[1].model_ymax[0] = ymax;
    double c0 = 1;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    int npx = atoi(argv[6]);
    int npy = atof(argv[7]);
    assert(npx > 1);
    assert(npy > 1);
    if (npx > 2) assert(npx % 2 != 0);
    if (npy > 2) assert(npy % 2 != 0);
    
    double theta = atof(argv[8]);
    assert(theta > -1e-6);
    theta *= 3.141592653589793 / 180.;
    
    double dz = atof(argv[9]);
    assert(dz > 0);
    
    double za = atof(argv[10]);
    double zb = atof(argv[11]);
    double zc = atof(argv[12]);
    int flag_3D = atoi(argv[13]);
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    // COUNT GRID ELEMENTS
    //++++++++++++++++++++++++++++++++++++++++++++++
    int nnodes = npx * npy;
    int nelems2d = 2*(npy-1)*(npx-1);
    int nelems1d = 2*(npx-1) + 2*(npy-1);
    //++++++++++++++++++++++++++++++++++++++++++++++
    
    FILE *fp,*fp_physics[nSupMods],*fp_params[nSupMods],*fp_design,*fp_cover,*fp_init[nSupMods];
    char *file_name_geo    = concat(file_name,".geo");
    char *file_name_design = concat(file_name,".design");
    char *file_name_cover  = concat(file_name,".cover");
    
    char *file_name_physics[nSupMods];
    char *file_name_params[nSupMods];
    char *file_name_init[nSupMods];
    char s1[50],s2[50], s3[50];
    for (i=0; i<nSupMods; i++) {
        sprintf(s1, "_sm%d.physics",i+1);
        file_name_physics[i] = concat(file_name,s1);
        fp_physics[i]=fopen(file_name_physics[i], "w");
        
        sprintf(s2, "_sm%d.params",i+1);
        file_name_params[i] = concat(file_name,s2);
        fp_params[i]=fopen(file_name_params[i], "w");
        for (j=0; j<sm[i].nMods; j++) {
            fprintf(fp_params[i],"MODEL %d %s\n",j+1,sm[i].model[j]);
        }
        fprintf(fp_params[i],"\n");
        fprintf(fp_params[i],"IP NTL 1e-6\n");
        fprintf(fp_params[i],"IP ITL 1e-6\n");
        fprintf(fp_params[i],"IP NIT 10\n");
        fprintf(fp_params[i],"IP MIT 100\n");
        
        sprintf(s3, "_sm%d.init",i+1);
        file_name_init[i] = concat(file_name,s3);
        fp_init[i]=fopen(file_name_init[i], "w");
        int itrns = 0
        for (j=0; j<sm[i].nMods; j++) {
            if (strcmp(sm[i].model[j], "SW2") == 0) {
                fprintf(fp_init[i],"DEPTH\n");
                for (k=0; k<nnodes; k++) {
                    fprintf(fp_init[i],"%lf\n",h0);
                }
                fprintf(fp_init[i],"ENDDS\n");
                fprintf(fp_init[i],"OLD_DEPTH\n");
                for (k=0; k<nnodes; k++) {
                    fprintf(fp_init[i],"%lf\n",h0);
                }
                fprintf(fp_init[i],"ENDDS\n");
                fprintf(fp_init[i],"OLDER_DEPTH\n");
                for (k=0; k<nnodes; k++) {
                    fprintf(fp_init[i],"%lf\n",h0);
                }
                fprintf(fp_init[i],"ENDDS\n");
                fprintf(fp_init[i],"DA-VELOCITY\n");
                for (k=0; k<nnodes; k++) {
                    fprintf(fp_init[i],"%lf %lf\n",u0,v0);
                }
                fprintf(fp_init[i],"ENDDS\n");
                fprintf(fp_init[i],"OLD_DA-VELOCITY\n");
                for (k=0; k<nnodes; k++) {
                    fprintf(fp_init[i],"%lf %lf\n",u0,v0);
                }
                fprintf(fp_init[i],"ENDDS\n");
                fprintf(fp_init[i],"OLDER_DA-VELOCITY\n");
                for (k=0; k<nnodes; k++) {
                    fprintf(fp_init[i],"%lf %lf\n",u0,v0);
                }
                fprintf(fp_init[i],"ENDDS\n");
            }
            if (strcmp(sm[i].model[j], "TRNS") == 0) {
                fprintf(fp_init[i],"CON %d\n");
                for (k=0; k<nnodes; k++) {
                    fprintf(fp_init[i],"%lf\n",c0);
                }
                fprintf(fp_init[i],"ENDDS\n");
                fprintf(fp_init[i],"OLD_CON %d\n");
                for (k=0; k<nnodes; k++) {
                    fprintf(fp_init[i],"%lf\n",c0);
                }
                fprintf(fp_init[i],"ENDDS\n");
                fprintf(fp_init[i],"OLDER_CON %d\n");
                for (k=0; k<nnodes; k++) {
                    fprintf(fp_init[i],"%lf\n",c0);
                }
                fprintf(fp_init[i],"ENDDS\n");
                itrns++;
            }
        }
    }
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Design File Write
    //++++++++++++++++++++++++++++++++++++++++++++++
    fp_design=fopen(file_name_design, "w");
    fprintf(fp_design,"#--------------------\n");
    fprintf(fp_design,"GRID %s\n",file_name_geo);
    fprintf(fp_design,"#--------------------\n");
    fprintf(fp_design,"COVERAGE %s\n",file_name_cover);
    fprintf(fp_design,"#--------------------\n");
    fprintf(fp_design,"# Time Controls\n");
    fprintf(fp_design,"#--------------------\n");
    fprintf(fp_design,"TC T0 0 0\n");
    fprintf(fp_design,"TC TF 100 0\n\n");
    fprintf(fp_design,"SERIES DT 1 2 0 0\n");
    fprintf(fp_design,"0 100\n");
    fprintf(fp_design,"100 100\n\n");
    fprintf(fp_design,"SERIES AWRITE 1 1 0 0\n");
    fprintf(fp_design,"0 86400 .1 0\n");
    fprintf(fp_design,"#--------------------\n");
    fprintf(fp_design,"# Coupled Model Setup\n");
    fprintf(fp_design,"#--------------------\n");
    for (i=0; i<nSupMods; i++) {
        fprintf(fp_design,"MONO %d %s_sm%d\n",i+1,file_name,i+1);
        if (i != 0) fprintf(fp_design,"COUPLE %d %d LAG\n",i,i+1);
    }
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    // Cover File Write
    //++++++++++++++++++++++++++++++++++++++++++++++
    // CJT - just assume all the same parameter material for now
    fp_cover=fopen(file_name_cover, "w");
    fprintf(fp_cover,"ELEMS_TRI %d\n",nelems2d);
    fprintf(fp_cover,"ELEMS_SEG %d\n",nelems1d);
    for (i=0; i<nelems2d; i++) {
        fprintf(fp_cover,"TRI %d\n",1);
    }
    for (i=0; i<nelems1d; i++) {
        fprintf(fp_cover,"SEG %d\n",1);
    }
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    // GEO File Write
    //++++++++++++++++++++++++++++++++++++++++++++++
    fp=fopen(file_name_geo, "w");
    fprintf(fp,"NODES %d\n",nnodes);
    fprintf(fp,"ELEMS_TRI %d\n",nelems2d);
    fprintf(fp,"ELEMS_SEG %d\n",nelems1d);
    for (i=0; i<nSupMods; i++) {
        fprintf(fp_physics[i],"ELEMS_TRI %d\n",nelems2d);
        fprintf(fp_physics[i],"ELEMS_SEG %d\n",nelems1d);
    }
    //++++++++++++++++++++++++++++++++++++++++++++++
    
    double dx = (xmax - xmin)/(double)(npx-1);
    double dy = (ymax - ymin)/(double)(npy-1);
    
    SVECT node[npx*npy];
    
    double x = xmin, xr, zr;
    double y = ymin, yr;
    k = 0;
    for (i=0; i<npx; i++) {
        for (j=0; j<npy; j++) {
            x = xmin + i*dx;
            y = ymax - j*dy; //ymin + j*dy;
            
            // totate 45 degrees
            xr = x*cos(theta) - y*sin(theta);
            yr = x*sin(theta) + y*cos(theta);
            zr = -(za + zb * x + zc * x * x);
            fprintf(fp,"ND %d %20.10f %20.10f %20.10f\n",k+1,xr,yr,zr);
            
            node[k].x = xr;
            node[k].y = yr;
            node[k].z = -(za + zb * x + zc * x * x);
            
            k++;
        }
    }
    assert(nnodes == k);
    
    SVECT cross;
    SVECT side1, side2;
    
    int inode=0, nd1, nd2, nd3, iMod, pmat;
    k=0;
    for (i=0; i<npx-1; i++) {
        for (igrp=0; igrp<npy-1; igrp++) {
            inode = i*npy + igrp;
            
            // make sure node numbering is counter clock-wise
            nd1 = inode; nd2 = inode+1; nd3 = inode+npy+1;
            side1.x =  node[nd2].x -  node[nd1].x; side1.y =  node[nd2].y -  node[nd1].y; side1.z =  node[nd2].z -  node[nd1].z;
            side2.x =  node[nd3].x -  node[nd1].x; side2.y =  node[nd3].y -  node[nd1].y; side2.z =  node[nd3].z -  node[nd1].z;
            cross = svect_cross(side1,side2);
            
            if (svect_dotp(cross,ref_vector) > 0) {
                //if (cross.z > 0) {
                nd1 = inode; nd2 = inode+npy+1; nd3 = inode+1;
                side1.x =  node[nd2].x -  node[nd1].x; side1.y =  node[nd2].y -  node[nd1].y; side1.z =  node[nd2].z -  node[nd1].z;
                side2.x =  node[nd3].x -  node[nd1].x; side2.y =  node[nd3].y -  node[nd1].y; side2.z =  node[nd3].z -  node[nd1].z;
                cross = svect_cross(side1,side2);
                assert(svect_dotp(cross,ref_vector) < 0);
                //assert(cross.z < 0);
            }
            //printf("node[0]: {%f,%f,%f}\n",node[nd1].x,node[nd1].y,node[nd1].z);
            //printf("node[1]: {%f,%f,%f}\n",node[nd2].x,node[nd2].y,node[nd2].z);
            //printf("node[2]: {%f,%f,%f}\n",node[nd3].x,node[nd3].y,node[nd3].z);
            //printf("nodes: %d %d %d \t cross: %20.10f %20.10f %20.10f \n",nd1+1,nd2+1,nd3+1,cross.x,cross.y,cross.z);
            
            // decide physics material
            for (iSupMod=0; iSupMod<nSupMods; iSupMod++) {
                pmat = -3;
                for (iMod=0; iMod<sm[iSupMod].nMods; iMod++) {
                    if ( (node[nd1].x >= sm[iSupMod].model_xmin[iMod] || node[nd2].x >= sm[iSupMod].model_xmin[iMod] || node[nd3].x >= sm[iSupMod].model_xmin[iMod]) &&
                         (node[nd1].x <= sm[iSupMod].model_xmax[iMod] || node[nd2].x <= sm[iSupMod].model_xmax[iMod] || node[nd3].x <= sm[iSupMod].model_xmax[iMod]) ) {
                        pmat = iMod;
                        fprintf(fp_physics[iSupMod],"TRI %d\n",pmat+1);
                        break;
                    }
                }
            }
            assert(pmat != -3);
            
            // original output
            fprintf(fp,"TRI %d %d %d %d -1\n",k+1,nd1+1,nd2+1,nd3+1); //,pmat+1);
            k++;
            
            
            // make sure node numbering is counter clock-wise
            nd1 = inode; nd2 = inode+npy+1; nd3 = inode+npy;
            side1.x =  node[nd2].x -  node[nd1].x; side1.y =  node[nd2].y -  node[nd1].y; side1.z =  node[nd2].z -  node[nd1].z;
            side2.x =  node[nd3].x -  node[nd1].x; side2.y =  node[nd3].y -  node[nd1].y; side2.z =  node[nd3].z -  node[nd1].z;
            cross = svect_cross(side1,side2);
            if (svect_dotp(cross,ref_vector) > 0) {
                //if (cross.z > 0) {
                nd1 = inode; nd2 = inode+npy; nd3 = inode+npy+1;
                side1.x =  node[nd2].x -  node[nd1].x; side1.y =  node[nd2].y -  node[nd1].y; side1.z =  node[nd2].z -  node[nd1].z;
                side2.x =  node[nd3].x -  node[nd1].x; side2.y =  node[nd3].y -  node[nd1].y; side2.z =  node[nd3].z -  node[nd1].z;
                cross = svect_cross(side1,side2);
                //assert(cross.z < 0);
                assert(svect_dotp(cross,ref_vector) < 0);
            }
            
            //printf("nodes: %d %d %d \t cross: %20.10f %20.10f %20.10f \n",nd1+1,nd2+1,nd3+1,cross.x,cross.y,cross.z);
            
            // decide physics material
            for (iSupMod=0; iSupMod<nSupMods; iSupMod++) {
                pmat = -3;
                for (iMod=0; iMod<sm[iSupMod].nMods; iMod++) {
                    if ( (node[nd1].x >= sm[iSupMod].model_xmin[iMod] || node[nd2].x >= sm[iSupMod].model_xmin[iMod] || node[nd3].x >= sm[iSupMod].model_xmin[iMod]) &&
                         (node[nd1].x <= sm[iSupMod].model_xmax[iMod] || node[nd2].x <= sm[iSupMod].model_xmax[iMod] || node[nd3].x <= sm[iSupMod].model_xmax[iMod]) ) {
                        pmat = iMod;
                        fprintf(fp_physics[iSupMod],"TRI %d\n",pmat+1);
                        break;
                    }
                }
            }
            assert(pmat != -3);
            
            fprintf(fp,"TRI %d %d %d %d -1\n",k+1,nd1+1,nd2+1,nd3+1); //,pmat+1);
            k++;
        }
    }
    assert(nelems2d == k);
    
    k=0;
    fprintf(fp,"! West Boundary\n");
    for (i=1; i<npy; i++) {
        for (iSupMod=0; iSupMod<nSupMods; iSupMod++) {
            fprintf(fp_physics[iSupMod],"SEG 1\n");
        }
        fprintf(fp,"SEG %d %d %d %d\n",k+1,i, i + 1, 5);
        k++;
    }
    
    fprintf(fp,"! South Boundary\n");
    for (i=1; i<npx; i++) {
        for (iSupMod=0; iSupMod<nSupMods; iSupMod++) {
            fprintf(fp_physics[iSupMod],"SEG 2\n");
        }
        fprintf(fp,"SEG %d %d %d %d\n",k+1,1 + npy * (i - 1), 1 + npy * i, 5);
        k++;
    }
    
    fprintf(fp,"! East Boundary\n");
    for (i=1; i<npy; i++) {
        for (iSupMod=0; iSupMod<nSupMods; iSupMod++) {
            fprintf(fp_physics[iSupMod],"SEG 3\n");
        }
        fprintf(fp,"SEG %d %d %d %d\n",k+1,1 + npy * (npx - 1) + (i - 1), 1 + npy * (npx - 1) + i, 5);
        k++;
    }
    
    fprintf(fp,"! North Boundary\n");
    for (i=1; i<npx; i++) {
        for (iSupMod=0; iSupMod<nSupMods; iSupMod++) {
            fprintf(fp_physics[iSupMod],"SEG 4\n");
        }
        fprintf(fp,"SEG %d %d %d %d\n",k+1,npy + npy * (i - 1), npy + npy * i, 5);
        k++;
    }
    assert(nelems1d == k);
    
    fclose(fp);
    
//    fp=fopen(file_name_bc, "w");
//    fprintf(fp,"# operational parameters\n");
//    fprintf(fp,"OP SW2\n");
//    fprintf(fp,"OP TRN 0\n");
//    fprintf(fp,"OP BLK 1\n");
//    fprintf(fp,"OP INC 40\n");
//    fprintf(fp,"OP PRE 1\n");
//    fprintf(fp,"\n");
//    fprintf(fp,"# solver parameters\n");
//    fprintf(fp,"IP NTL 1e-6\n");
//    fprintf(fp,"IP ITL 1e-6\n");
//    fprintf(fp,"IP NIT 10\n");
//    fprintf(fp,"IP MIT 100\n");
//    fprintf(fp,"\n");
//    fprintf(fp,"# 2d element material string\n");
//    fprintf(fp,"MTS 1 1\n");
//    fprintf(fp,"\n");
//    fprintf(fp,"# material properties\n");
//    fprintf(fp,"MP ML 1 0\n");
//    fprintf(fp,"MP SRT 1 100\n");
//    fprintf(fp,"MP EVS 1 0.0 0.0 0.0\n");
//    fprintf(fp,"MP MUC 1.0\n");
//    fprintf(fp,"MP MU 1e-6\n");
//    fprintf(fp,"MP RHO 1000\n");
//    fprintf(fp,"MP G 9.8\n");
//    fprintf(fp,"\n");
//    fprintf(fp,"# friction coefficients\n");
//    fprintf(fp,"FR MNG 1 0.0\n");
//    fprintf(fp,"\n");
//    fprintf(fp,"#TIME-SERIES --------------------\n\n");
//    fprintf(fp,"# Output time-series \n");
//    fprintf(fp,"SERIES AWRITE 1 1 0 0\n");
//    fprintf(fp,"%-20.10f %-20.10f %-20.10f 0\n",0.,1000.,100.);
//    fprintf(fp,"\n");
//    fprintf(fp,"# Time-step time-series\n");
//    fprintf(fp,"SERIES DT 2 2 0 0\n");
//    fprintf(fp,"%-20.10f %-20.10f\n",0.,100.);
//    fprintf(fp,"%-20.10f %-20.10f\n",1000.,100.);
//    fprintf(fp,"\n");
//    fprintf(fp,"# No-flow time-series\n");
//    fprintf(fp,"SERIES BC 3 2 0 0\n");
//    fprintf(fp,"%-20.10f 0.0\n",0.);
//    fprintf(fp,"%-20.10f 0.0\n",1000.);
//    fprintf(fp,"\n");
//    fprintf(fp,"TC T0 %-20.10f\n",0.);
//    fprintf(fp,"TC TF %-20.10f\n",1000.);
//    fprintf(fp,"\n");
//    fprintf(fp,"#STRINGS ------------------------\n");
//    fprintf(fp,"! West\n");
//    fprintf(fp,"NB VEL 2 3\n");

//    fprintf(fp,"\nEND\n");
//    fclose(fp);
//    printf("bc file written\n");
//
//
//    fp=fopen(file_name_hot, "w");
//    fprintf(fp,"DATASET\n");
//    fprintf(fp,"OBJTYPE \"mesh2d\"\n");
//    fprintf(fp,"BEGSCL\n");
//    fprintf(fp,"ND %d\n",npoints);
//    fprintf(fp,"NC %d\n",nelemsTotal);
//    fprintf(fp,"NAME ioh\n");
//    fprintf(fp,"TS 0 0\n");
//    for (i=0; i<npx; i++) {
//        for (j=0; j<npy; j++) {
//            x = xmin + i*dx;
//            y = ymax - j*dy; //ymin + j*dy;
//
//            // totate 45 degrees
//            xr = x*cos(theta) - y*sin(theta);
//            yr = x*sin(theta) + y*cos(theta);
//            zr = -(za + zb * x + zc * x * x);
//            fprintf(fp,"%-20.10f\n",-zr);
//        }
//    }
//    fprintf(fp,"ENDDS\n");
//    fclose(fp);
//    printf("hotstart written\n");
//
//
//    if (flag_3D) {
//        fp=fopen(file_name_bin, "w");
//        fprintf(fp,"BIN 1 1 0. %f \n",dz);
//        fprintf(fp,"BIN 1 2 -1000 %f \n",dz);
//        fclose(fp);
//        printf("extrusion bin file written\n");
//
//        fp=fopen(file_name_node, "w");
//        for (i=0; i< npoints; i++) {
//            fprintf(fp,"ND %d 1 \n",i+1);
//        }
//        fclose(fp);
//        printf("extrusion node file written\n");
//
//
//        free(file_name_2dm);
//        free(file_name_bc);
//        free(file_name_node);
//        free(file_name_hot);
//        free(file_name_bin);
//
//
//        // NOW CALL 3D EXTRUSION
//        printf("extruding files now ... \n");
//        int flag_bins = 1;
//        int flag_min_layer = 0;
//        int minimum_layers = 1;
//        int mesh_type = TETRAHEDRON; //MIXED_ELEMENT_MESH;
//        extrudeAdH(flag_bins,flag_min_layer,minimum_layers,file_name,mesh_type);
//
//    }
    
    return 0;
}


