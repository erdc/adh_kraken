#ifndef H_XDMF_
#define H_XDMF_

void write_xdmf_header(FILE *xmf);
void write_xdmf_tail(FILE *xmf);
void init_hdf5_file(MPI_Comm comm, MPI_Info info, int mpi_rank, int mpi_size, char *fbase);
void write_hdf5_mesh(MPI_Comm comm, MPI_Info info, int mpi_rank, int mpi_size, char *fbase, int nt);
void write_hdf5_data(MPI_Comm comm, MPI_Info info, int mpi_rank, int mpi_size, char *fbase, int nt, float t);
void xdmf_init(FILE *xmf, char *fbase);
void xdmf_write_dataset(FILE *xmf, char *fbase, int nt, int mesh_no, float t);
void xdmf_finalize(FILE *xmf, char *fbase);

#endif
