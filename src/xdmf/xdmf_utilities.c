/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  This file contains a collection of functions for writing xmf files                                                                                    */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
#define RANK 2
#define SPATIAL_DIM 3 // Always be 3
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Initializes and XMF file and writes the header
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] fbase   (char *)  the xmf filename
 * @param[in] dn          (char *)  the xmf domain name
 * @param[in] gn          (char *)  the xmf grid name
 * \return xmf (FILE \*) a file pointer to an XMF file
 *
 * \note CJT - only works for serial at the moment
 * \note Tutorial at: https://docs.hdfgroup.org/hdf5/v1_14/v1_14_4/_ex_a_p_i.html
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
FILE *xmf_init(char *fbase, char *dn, char *gn) {
#ifdef _HDF5
    char fn_xmf[50];
    strcpy(fn_xmf,fbase);
    strcat(fn_xmf, ".xmf");
    FILE *xmf = fopen(fn_xmf, "w");
    xmf_header(xmf,dn,gn);
    return xmf;
#else
    return NULL;
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Finalized an XMF file and writes end matter
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] xmf   (FILE *)  the xmf file
 *
 * \note CJT - only works for serial at the moment
 * \note Tutorial at: https://docs.hdfgroup.org/hdf5/v1_14/v1_14_4/_ex_a_p_i.html
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void xmf_finalize(FILE *xmf) {
#ifdef _HDF5
    xmf_tail(xmf);
    fclose(xmf);
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes an XMF header
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] xmf   (FILE *)  the xmf file
 * @param[in] dn          (char *)  the xmf domain name
 * @param[in] gn          (char *)  the xmf grid name
 *
 * \note CJT - only works for serial at the moment
 * \note Tutorial at: https://docs.hdfgroup.org/hdf5/v1_14/v1_14_4/_ex_a_p_i.html
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void xmf_header(FILE *xmf, char *dn, char *gn){
#ifdef _HDF5
#ifdef _MESSG
    if(myid==0)
#endif
    {
        fprintf(xmf, "<Xdmf Version=\"3.0\">\n");
        fprintf(xmf, "\t <Domain Name=\"%s\">\n",dn);
        fprintf(xmf, "\t\t<Grid Name=\"%s\" GridType=\"Collection\" CollectionType=\"Temporal\">\n",gn);
    }
#endif
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes XMF file end matter
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] xmf   (FILE *)  the xmf file
 *
 * \note CJT - only works for serial at the moment
 * \note Tutorial at: https://docs.hdfgroup.org/hdf5/v1_14/v1_14_4/_ex_a_p_i.html
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void xmf_tail(FILE *xmf){
#ifdef _HDF5
#ifdef _MESSG
    if(myid==0)
#endif
    {
        fprintf(xmf, "\t\t</Grid>\n");
        fprintf(xmf, "\t</Domain>\n");
        fprintf(xmf, "</Xdmf>\n");
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes XMF time-step header
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] myid   (int)  my processor ID
 * @param[in] xmf     (FILE *)  the xmf file
 * @param[in] fbase (char *)  the xmf grid name
 * @param[in] time (double) the current model time
 *
 * \note CJT - only works for serial at the moment
 * \note Tutorial at: https://docs.hdfgroup.org/hdf5/v1_14/v1_14_4/_ex_a_p_i.html
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void xmf_write_ts_header(int myid, FILE *xmf, char *fbase, double time) {
#ifdef _HDF5
#ifdef _MESSG
    if(myid==0)
#endif
    {
        fprintf(xmf, "\t\t\t<Grid Name=\"%s\">\n",fbase);
        fprintf(xmf, "\t\t\t<Time Value=\"%lf\" />\n",time);
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes XMF time-step tail
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] xmf     (FILE *)  the xmf file
 *
 * \note CJT - only works for serial at the moment
 * \note Tutorial at: https://docs.hdfgroup.org/hdf5/v1_14/v1_14_4/_ex_a_p_i.html
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void xmf_write_ts_tail(FILE *xmf) {
#ifdef _HDF5
#ifdef _MESSG
    if(myid==0)
#endif
    {
        fprintf(xmf, "\t\t\t</Grid>\n");
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes XMF grid geometry
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] g              (SGRID *)  and AdH grid
 * @param[in] xmf          (FILE *)  the xmf file
 * @param[in] fbase      (char *)  the xmf filename
 * @param[in] time        (double) the current model time
 * @param[in] mesh_no (int) the adapted mesh grid number.  This is 0 if not printing adapted meshes.
 *
 * \note CJT - only works for serial at the moment
 * \note Tutorial at: https://docs.hdfgroup.org/hdf5/v1_14/v1_14_4/_ex_a_p_i.html
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void xmf_write_geometry(SGRID *g, FILE *xmf, char *fbase, int mesh_no) {
#ifdef _HDF5
#ifdef _MESSG
    if(g->smpi->myid==0)
#endif
    {
        fprintf(xmf, "\t\t\t\t<Geometry GeometryType=\"XYZ\">\n");
        fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d %d\" Format=\"HDF\">\n", g->macro_nnodes,SPATIAL_DIM);
        fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Mesh/XY/%d\n",fbase,mesh_no);
        fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
        fprintf(xmf, "\t\t\t\t</Geometry>\n");
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes XMF grid topologu
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] g              (SGRID *)  and AdH grid
 * @param[in] xmf          (FILE *)  the xmf file
 * @param[in] fbase      (char *)  the xmf filename
 * @param[in] time        (double) the current model time
 * @param[in] mesh_no (int) the adapted mesh grid number.  This is 0 if not printing adapted meshes.
 *
 * \note CJT - only works for serial at the moment
 * \note Tutorial at: https://docs.hdfgroup.org/hdf5/v1_14/v1_14_4/_ex_a_p_i.html
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void xmf_write_topology(SGRID *g, FILE *xmf, char *fbase, int mesh_no) {
#ifdef _HDF5
#ifdef _MESSG
    if(g->smpi->myid==0)
#endif
    {
        int macro_nelems = g->macro_nelems1d + g->macro_nelems2d + g->macro_nelems3d;
        int nentry = g->macro_nSegs * 3 +  g->macro_nQuads * 5 + g->macro_nTris * 4 + g->macro_nTets * 5 + g->macro_nPrisms * 7; // cjt -- this should be fine
        fprintf(xmf, "\t\t\t\t<Topology TopologyType=\"Mixed\" NumberOfElements=\"%d\">\n",macro_nelems);
        
        nentry = g->macro_nQuads * 5 + g->macro_nTris * 4 + g->macro_nTets * 5 + g->macro_nPrisms * 7;  // cjt -- not sure why I can't write out the segments.
        
        fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"HDF\">\n",nentry);
        fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Mesh/Elements/%d\n",fbase,mesh_no);
        fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
        fprintf(xmf, "\t\t\t\t</Topology>\n");
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Writes XMF data attribute
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] g              (SGRID *)  and AdH grid
 * @param[in] xmf          (FILE *)  the xmf file
 * @param[in] fbase      (char *)  the xmf filename
 * @param[in] name        (char *)  the variable name
 * @param[in] ndim        (int) the dimensionality of the data
 * @param[in] nt             (int) the AdH time-step
 * @param[in] isNodal (bool) true - node-based data, false - element-based data
 *
 * \note CJT - only works for serial at the moment
 * \note Tutorial at: https://docs.hdfgroup.org/hdf5/v1_14/v1_14_4/_ex_a_p_i.html
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void xmf_write_ts_attribute(SGRID *g, FILE *xmf, char *fbase, char *name, int ndim, int nt, bool isNodal) {
#ifdef _HDF5
#ifdef _MESSG
    if(g->smpi->myid==0)
#endif
    {
        bool isScalar = false;  if (ndim == 1) isScalar = true;
        if (isNodal) {
            if (isScalar) {
                fprintf(xmf, "\t\t\t\t<Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",name);
                fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", g->macro_nnodes);
                fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/%s/%d\n",fbase,name,nt);
            } else {
                fprintf(xmf, "\t\t\t\t<Attribute Name=\"%s\" AttributeType=\"Vector\" Center=\"Node\">\n",name);
                fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d %d\" Format=\"HDF\">\n", g->macro_nnodes,ndim);
                fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/%s/%d\n",fbase,name,nt);
            }
        } else {
            if (isScalar) {
                fprintf(xmf, "\t\t\t\t<Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",name);
                fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", g->macro_nelems1d + g->macro_nelems2d + g->macro_nelems3d);
                fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/%s/%d\n",fbase,name,nt);
            } else {
                fprintf(xmf, "\t\t\t\t<Attribute Name=\"%s\" AttributeType=\"Vector\" Center=\"Cell\">\n",name);
                fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d %d\" Format=\"HDF\">\n", g->macro_nelems1d + g->macro_nelems2d + g->macro_nelems3d,ndim);
                fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/%s/%d\n",fbase,name,nt);
            }
        }
        fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
        fprintf(xmf, "\t\t\t\t</Attribute>\n");
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
