#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <hdf5.h>
#include <mpi.h>
#include <assert.h>

void xdmf_init(FILE *xmf, char *fbase){
    char fname[50];
    strcpy(fname, fbase);
    strcat(fname, ".xmf");
    xmf = fopen(fname, "w");
    //breaks stuff?
    //fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    //fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"3.0\">\n");
    //Domain information
    fprintf(xmf, "\t <Domain Name=\"Adh Sim\">\n");
    //Time based grid start
    fprintf(xmf, "\t\t<Grid Name=\"MeshTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");
    fclose(xmf);
}
