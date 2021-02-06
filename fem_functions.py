import numpy as np
from os import linesep as endl

from data import *

def Feeldof(nd, nnel, ddln):
    edof = nnel*ddln
    index = np.zeros(edof, dtype=np.int32)
    k=0
    for i in range(0,nnel):
        start = nd[i]*ddln
        for j in range(0,ddln):
            index[k] = start+j
            k=k+1
    return(index)

def Assembl(KG,ke,index):
    edof = len(index)
    for i in range(0,edof):
        ii=index[i]
        for j in range(0,edof):
            jj=index[j]
            KG[ii,jj]=KG[ii,jj]+ke[i,j]
    return(KG)

def AssemblF(FG,fe,index):
    edof = len(index)
    for i in range(0,edof):
        ii=index[i]
        FG[ii]=FG[ii]+fe[i]
    return(FG)

def Interpolc(ksi):
    # Cette fonction évalue numériquement la valeur des fonctions N aux
    # coordonnées demandées en Ksi:
    N = np.array([(1-ksi)/2., (1+ksi)/2.]);
    return(N)

def Interpol(ksi,eta):
    # Cette fonction évalue numériquement la valeur des fonctions N aux
    # coordonnées demandées en Ksi et Eta:
    N = np.array([(1-ksi-eta+ksi*eta)/4.,(1+ksi-eta-ksi*eta)/4.,(1+ksi+eta+ksi*eta)/4.,(1-ksi+eta-ksi*eta)/4.], dtype=float_type)
    return(N)

def BKsiEta(ksi,eta):
    Bke = np.array([[(eta-1)/4., (1-eta)/4., (1+eta)/4., -(1+eta)/4.],[(ksi-1)/4., -(1+ksi)/4., (1+ksi)/4., (1-ksi)/4.]], dtype=float_type)
    return(Bke)

def Gauss(npg):
    if(npg == 2):
        xloc = np.array([1/sqrt(3.), -1./sqrt(3)], dtype=float_type)
        W = np.array([1., 1.], dtype=float_type);
    if(npg == 4):
        # Méthode produit à 2 X 2 points
        # 4 points intègrent exactement jusqu'à l'ordre 3
        W = np.array([1, 1, 1, 1], dtype=float_type)
        c = 1./sqrt(3.)
        xloc = np.array([[c, c], [-c, c], [-c, -c], [c, -c]])
    if(npg == 9):
        #Points et poids de Gauss pour un carré
        w1 = np.array([8/9., 5/9., 5/9.], dtype=float_type)
        x1 = np.array([0, -sqrt(3/5.), sqrt(3/5.)])
        # Méthode à 9 points intégrant jusqu'à l'ordre m=5 dans chaque direction
        W = np.zeros(9, dtype=float_type)
        xloc = np.zeros((9,2), dtype=float_type)
        k=0;
        for i in range(0,3):
            for j in range(0,3):
                W[k]=w1[i]*w1[j];
                xloc[k,0]=x1[i];
                xloc[k,1]=x1[j];
                k=k+1;
    return(xloc,W)
