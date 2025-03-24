#ifndef H_XDMF_
#define H_XDMF_

FILE *xmf_init(char *fbase, char *dn, char *gn, int myid);
void xmf_finalize(FILE *xmf, int myid);
void xmf_header(FILE *xmf, char *dn, char *gn, int myid);
void xmf_tail(FILE *xmf, int myid);
void xmf_write_ts_header(int myid, FILE *xmf, char *fbase, double time);
void xmf_write_ts_tail(FILE *xmf, int myid);
void xmf_write_geometry(SGRID *g, FILE *xmf, char *fbase, int mesh_no);
void xmf_write_topology(SGRID *g, FILE *xmf, char *fbase, int mesh_no);
void xmf_write_ts_attribute(SGRID *g, FILE *xmf, char *fbase, char *name, int ndim, int nt, bool isNodal);



void hdf5_init(SGRID *g, char *fbase);
void hdf5_write_data(SGRID *g, char *fbase, char *name, double *data, int ndim, bool isNodal, int nt);
void hdf5_write_grid(SGRID *g, char *fbase);
// so does data write


//void write_xdmf_header(FILE *xmf);
//void write_xdmf_tail(FILE *xmf);
//FILE *init_hdf5_file(SGRID *g, char *fbase);
//void xdmf_init(FILE *xmf, char *fbase);
//void xdmf_finalize(FILE *xmf, char *fbase);
//void write_hdf5_data(SGRID *g, FILE *xmf, char *fbase, char *name, double *data, int ndim, bool isNodal, int nt);

#endif
