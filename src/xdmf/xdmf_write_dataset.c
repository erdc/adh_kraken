#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <hdf5.h>
#include <mpi.h>
#include <assert.h>
#define NumElements 4
#define NumNodes 8
#define NEntry 18
#define RANK 2
//RANK is dimension of data set, 2 is a matrix, 1 is a vector

void xdmf_write_dataset(FILE *xmf, char *fbase, int nt, int mesh_no, float t){
    char fname[50];
    strcpy(fname, fbase);
    strcat(fname, ".xmf");
    xmf = fopen(fname, "a");

    //Grid start
        fprintf(xmf, "\t\t\t<Grid Name=\"2D Unstructured Mesh\">\n");
            fprintf(xmf, "\t\t\t<Time Value=\"%f\" />\n",t);
            //Geometry start, this is nodes
            fprintf(xmf, "\t\t\t\t<Geometry GeometryType=\"XY\">\n");
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d %d\" Format=\"HDF\">\n", NumNodes,2);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Mesh/XY/%d\n",fbase,mesh_no);
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Geometry>\n");
            //Geometry finish, this is nodes

            //Topology start
            fprintf(xmf, "\t\t\t\t<Topology TopologyType=\"Mixed\" NumberOfElements=\"%d\">\n", NumElements);
            //Data Item is the mixed connectivity
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", NEntry);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Mesh/Elements/%d\n",fbase,mesh_no);    
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Topology>\n");
            //Topology finish, this contains connectivities
            
            
            //Any nodal data we can also write
            fprintf(xmf, "\t\t\t\t<Attribute Name=\"NodalData\" AttributeType=\"Scalar\" Center=\"Node\">\n");
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", NumNodes);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/NodalScalar/%d\n",fbase,nt);
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Attribute>\n");
            //Nodal data finish
            

            //Should also be way to write elemental data here too just change Center=\"Node" to Center=\"Cell"
            fprintf(xmf, "\t\t\t\t<Attribute Name=\"ElementalData\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", NumElements);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/ElementalScalar/%d\n",fbase,nt);
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Attribute>\n");
            
            
            //Should also be way to write elemental data here too just change Center=\"Node" to Center=\"Cell"
            fprintf(xmf, "\t\t\t\t<Attribute Name=\"NodalVector\" AttributeType=\"Vector\" Center=\"Node\">\n");
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d %d\" Format=\"HDF\">\n", NumNodes,RANK);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/NodalVector/%d\n",fbase,nt);
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Attribute>\n");

            //Should also be way to write elemental data here too just change Center=\"Node" to Center=\"Cell"
            fprintf(xmf, "\t\t\t\t<Attribute Name=\"ElementalVector\" AttributeType=\"Vector\" Center=\"Cell\">\n");
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d %d\" Format=\"HDF\">\n", NumElements,RANK);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/ElementalVector/%d\n",fbase,nt);
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Attribute>\n");
            
            
            
        
        fprintf(xmf, "\t\t\t</Grid>\n");
        fclose(xmf);
}
