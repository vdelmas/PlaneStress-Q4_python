import numpy as np
from os import linesep as endl

from data import *

class mesh_cl():
    def __init__(self,nodesfile,elemsfile,nnel,contourfile,nnelc):
        self.nnel = nnel # Nombre de noeud par élément intérieur
        self.nnelc = nnelc # Nombre de noeuds par element de contour

        #Lecture du fichier de noeud
        self.nodes = np.loadtxt(nodesfile, dtype=float_type)[1:,1:] # Enlève la première ligne & collone
        self.nnodes = len(self.nodes)

        #Lecture du fichier élements interne
        self.elems = np.loadtxt(elemsfile, dtype=np.int32)[1:,11:]-1
        self.nelems = len(self.elems)
        assert(len(self.elems[0])==nnel)

        #Lecture du fichier élements contour
        self.elems_contour = np.loadtxt(contourfile, dtype=np.int32)[1:,1:]-1
        self.nelems_contour = len(self.elems_contour)
        assert(len(self.elems_contour[0])==nnelc)


def export_mesh(mesh):        
    filename = "mesh.vtk"
    with open(filename,"w") as f:
        f.write("# vtk DataFile Version 1.0" + endl)
        f.write("PlaneStress Mesh SYS806" + endl)
        f.write("ASCII" + endl)
        f.write("DATASET UNSTRUCTURED_GRID" + endl)
        f.write("POINTS " +  str(mesh.nnodes) + " float" + endl)
        for i in range(0,mesh.nnodes):
            f.write(str(mesh.nodes[i][0]) + " " + str(mesh.nodes[i][1]) + " " + str(mesh.nodes[i][2]) + endl)
        f.write(endl + "CELLS " + str(mesh.nelems+mesh.nelems_contour) + " " + str(mesh.nelems*(mesh.nnel+1)+mesh.nelems_contour*(mesh.nnelc+1)) + endl)
        for i in range(0,mesh.nelems):
            f.write(str(mesh.nnel) + " " ) 
            for j in range(mesh.nnel):
                f.write(str(mesh.elems[i][j]) + " ")
            f.write(endl)
        for i in range(0,mesh.nelems_contour):
            f.write(str(mesh.nnelc) + " " ) 
            for j in range(mesh.nnelc):
                f.write(str(mesh.elems_contour[i][j]) + " ")
            f.write(endl)

        f.write(endl + "CELL_TYPES " + str(mesh.nelems+mesh.nelems_contour) + endl)
        if(mesh.nnel == 3):
            cellType_elem=5; #Triangle
        elif(mesh.nnel == 4):
            cellType_elem=9; #Quadrangle
        else:
            raise Exception("ELEM CELL TYPE ERROR");
        if(mesh.nnelc == 2):
            cellType_elem_contour=3; #Line
        else:
            raise Exception("ELEM CONTOUR CELL TYPE ERROR");

def export_sol_vtk(mesh, u, v, sigmax, sigmay, sigmaxy):        
    filename = "mesh_sol.vtk"
    with open(filename,"w") as f:
        f.write("# vtk DataFile Version 1.0" + endl)
        f.write("PlaneStress Mesh SYS807" + endl)
        f.write("ASCII" + endl)
        f.write("DATASET UNSTRUCTURED_GRID" + endl)
        f.write("POINTS " +  str(mesh.nnodes) + " float" + endl)
        for i in range(0,mesh.nnodes):
            f.write(str(mesh.nodes[i][0]) + " " + str(mesh.nodes[i][1]) + " " + str(mesh.nodes[i][2]) + endl)

        f.write(endl + "CELLS " + str(mesh.nelems+mesh.nelems_contour) + " " + str(mesh.nelems*(mesh.nnel+1)+mesh.nelems_contour*(mesh.nnelc+1)) + endl)
        for i in range(0,mesh.nelems):
            f.write(str(mesh.nnel) + " " ) 
            for j in range(mesh.nnel):
                f.write(str(mesh.elems[i][j]) + " ")
            f.write(endl)
        for i in range(0,mesh.nelems_contour):
            f.write(str(mesh.nnelc) + " " ) 
            for j in range(mesh.nnelc):
                f.write(str(mesh.elems_contour[i][j]) + " ")
            f.write(endl)

        f.write(endl + "CELL_TYPES " + str(mesh.nelems+mesh.nelems_contour) + endl)
        if(mesh.nnel == 3):
            cellType_elem=5; #Triangle
        elif(mesh.nnel == 4):
            cellType_elem=9; #Quadrangle
        else:
            raise Exception("ELEM CELL TYPE ERROR");
        if(mesh.nnelc == 2):
            cellType_elem_contour=3; #Line
        else:
            raise Exception("ELEM CONTOUR CELL TYPE ERROR");

        for i in range(0,mesh.nelems):
            f.write(str(cellType_elem) + endl) 
        for i in range(0,mesh.nelems_contour):
            f.write(str(cellType_elem_contour) + endl) 

        f.write(endl + "POINT_DATA " + str(mesh.nnodes) + endl)
        f.write("Vectors Deplacement_u_v_(mm) float" + endl)
        for i in range(0,mesh.nnodes):
            f.write(str(u[i])+" "+str(v[i])+" 0.0" + endl)
        f.write("Scalars SigmaX(MPa) float" + endl)
        f.write("LOOKUP_TABLE default" + endl)
        for i in range(0,mesh.nnodes):
            f.write(str(sigmax[i]) + endl)
        f.write("Scalars SigmaY(MPa) float" + endl)
        f.write("LOOKUP_TABLE default" + endl)
        for i in range(0,mesh.nnodes):
            f.write(str(sigmay[i]) + endl)
        f.write("Scalars Sigmaxy(MPa) float" + endl)
        f.write("LOOKUP_TABLE default" + endl)
        for i in range(0,mesh.nnodes):
            f.write(str(sigmaxy[i]) + endl)

