#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <hdf5.h>
#include <mpi.h>
#include <assert.h>

void xdmf_finalize(FILE *xmf, char *fbase){
    char fname[50];
    strcpy(fname,fbase);
    strcat(fname, ".xmf");
    xmf = fopen(fname, "a");


     //Grid Finish
    fprintf(xmf, "\t\t</Grid>\n");
    fprintf(xmf, "\t</Domain>\n");
    //Domain finished
    fprintf(xmf, "</Xdmf>\n");
    //XDMF File finished
    
    fclose(xmf);

}
