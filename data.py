from math import *
import numpy as np

#simple vs double precision
float_type = np.float64

E = 2*10**11     # Module de Young en Pa
v = 0.3         # Coefficient de Poisson
Hypothesis = 1  # Plane stress
rho = 7850      # Masse volumique en Kg/m3
b = 0.02        # épaisseur en m selon z

# Facteur d'amplification pour la visualisation
ampli = 100;

# Charge répartie 
q = 10**6;    # N/m répartis uniformément de x=0 à x=L

# Forces volumiques (poids, inertie, force magnétique) (m/s**2)
fx = 0;       
fy = 0;

# Forces surfaciques (N/m**2)'
L = sqrt( 0.6**2+ 1.5**2);
sin = 1.5/L;
cos =  0.6/L;
fsx = -(q/b) *(sin)
fsy =  (q/b) *(cos)

# Lecture du nombre de ddl par noeud
ddln = 2;
# Lecture du nombre de ddl par noeud sur les contours
ddlnc = 2;
# Lecture du nombre de noeud par élément interne
nnel = 4;
# Nombre de noeud de l'élément de contour
nnelc = 2;
# Calcul du nombre de DDL par élément interne (8 dans le cas du Q4)
Nui = ddln*nnel;
# Calcul du nombre de DDL par élément de contour (4 pour barre linéaire)
Nuic = ddlnc*nnelc;

# Matrice de propriétés physiques
if(Hypothesis  == 1):    # plane stress
    facteur = E/(1-v**2);
    D = facteur*np.array([[1 , v, 0],[v , 1, 0],[0 , 0, (1-v)/2]], dtype=float_type);
else:
    facteur = E/((1-2*v)*(1+v));  # plane strain
    # D = facteur*np.array([[1-v, v, 0],[v , 1-v, 0],[0, 0, (1-2*v)/2.]]); ??
    mu = E/(2+2*v);
    D1 = np.array([[2*mu, 0, 0],[0, 2*mu, 0],[0, 0, mu]], dtype=float_type);
    lambd = facteur*v;
    D2 = facteur*np.array([[1, 1, 0],[1, 1, 0],[0, 0, 0]], dtype=float_type);
    D = D1+D2;
