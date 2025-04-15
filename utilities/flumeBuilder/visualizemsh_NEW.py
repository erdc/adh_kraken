import h5py as h5
import sys
import numpy as np
import meshio
import matplotlib.pyplot as plt
import matplotlib.tri as tri
def read_adh_mesh(fname):
    f = open(fname,"r")
    tets = []
    prisms = []
    tris = []
    quads = []
    nds = []
    for line in f:
        linelist = line.split()
        if linelist[0] == "TET":
            tet_no = int(linelist[1])
            tet_nodes = list(map(int,linelist[2:6]))
            tets.append(tet_nodes)
        elif linelist[0] =='PRISM':
            prism_no = int(linelist[1])
            prism_nodes = list(map(int,linelist[2:8]))
            prisms.append(prism_nodes)
        elif linelist[0] == "TRI":
            tri_no = int(linelist[1])
            tri_nodes = list(map(int,linelist[2:5]))
            tris.append(tri_nodes)
        elif linelist[0] == "QUAD":
            quad_no = int(linelist[1])
            quad_nodes = list(map(int,linelist[2:6]))
            quads.append(quad_nodes)
        elif linelist[0] == "ND":
            node_no = int(linelist[1])
            nd_coords = list(map(float,linelist[2:5]))
            nds.append(nd_coords)
            
    f.close()
    return tets,prisms,tris,quads,nds


tets,prisms,tris,quads,nds = read_adh_mesh(sys.argv[1])

#adjust numbers to be starting at 0
tets = np.array(tets,dtype=np.int32)
tris = np.array(tris,dtype=np.int32)
tets  = tets-1
tris = tris-1



cells = [
    ("triangle", tris)#,
#    ("tetra", tets)
]

mesh = meshio.Mesh(
    nds,
    cells
    # Optionally provide extra data on points, cells, etc.
    #point_data={"T": [0.3, -1.2, 0.5, 0.7, 0.0, -3.0]},
    # Each item in cell data must match the cells array
    #cell_data={"a": [[0.1, 0.2], [0.4]]},
)


mesh.write(
    "newgeo_test.vtk",  # str, os.PathLike, or buffer/open file
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)

mesh.write("newgeo_test.xmf")

#also plot on matplotlib
nds = np.array(nds)
tris = tri.Triangulation(nds[:,0],nds[:,1],triangles=tris)
plt.triplot(tris)
plt.savefig("Mesh.png")
