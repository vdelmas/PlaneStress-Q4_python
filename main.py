import numpy as np

from data import *          # Données du problème
from fem_functions import *  # Classes de la méthode fem
from input_output import *  # Fonctions liées au maillage et à l'export de la solution

#Lecture du maillage via la classe mesh_cl
mesh = mesh_cl("nodes.txt","elements.txt",4,"contour.txt",2)
export_mesh(mesh)

nnt = mesh.nnodes
nelt = mesh.nelems
nelc = mesh.nelems_contour
ndlt = nnt*ddln;

#Lecture des conditions aux limites 
CLx_0 = np.loadtxt("DirichletX.txt",dtype=np.int32)[1:,0]-1 #Correction de la numérotation pour commencer à 0
CLx_1 = np.loadtxt("DirichletX.txt",dtype=np.int32)[1:,1]
nclx = len(CLx_0)

CLy_0 = np.loadtxt("DirichletY.txt",dtype=np.int32)[1:,0]-1 #Correction de la numérotation pour commencer à 0
CLy_1 = np.loadtxt("DirichletY.txt",dtype=np.int32)[1:,1]
ncly = len(CLy_0)

#Initialisation des matrices 
#peut être modifié pour utiliser le format scipy.sparse
KG= np.zeros((ndlt,ndlt), dtype=float_type)
MG = np.zeros((nnt,nnt), dtype=float_type)
KN = np.zeros((ndlt,ndlt), dtype=float_type)

FG = np.zeros(ndlt, dtype=float_type);
RG = np.zeros(ndlt, dtype=float_type);

#Définition de l'intégration de Gauss
npg = 4; 
(xloc,W) = Gauss(npg);

#Boucle sur les éléments
for ie in range(0, nelt):
    # Extraire les infos propres à cet élément (coord., connect.,propriétés. ...)
    # Boucle sur les noeuds de l'élément en question 
    nd = np.zeros(nnel, dtype=np.int32)
    coor = np.zeros((nnel,2), dtype=float_type) # Removing z coordinate
    for i in range(0,nnel):
        nd[i] = mesh.elems[ie,i];
        coor[i,0] = mesh.nodes[nd[i],0];  # Coordonnées des noeuds de l'élément
        coor[i,1] = mesh.nodes[nd[i],1];  # Coordonnées des noeuds de l'élément

    # Initialisation des matrices et vecteurs élémentaires
    ke = np.zeros((Nui,Nui), dtype=float_type);
    me = np.zeros((nnel,nnel), dtype=float_type);
    fe = np.zeros(Nui, dtype=float_type);

    # Boucle sur les points d'intégration
    for i in range(0,npg):
        # Extraction des coordonnees des points d'intégration
        k = xloc[i,0];  # Ksi du point d'intégration sur le carré en cours
        n = xloc[i,1];  # Eta du point d'intégration sur le carré en cours

        # Évaluation numérique de la matrice Bksieta pour le point et l'élément
        Bke = BKsiEta(k,n);
        # Évaluation numérique des fonctions N pour le point et l'élément
        Nnum = Interpol(k,n);
            
        # Calcul numérique de la matrice Jacobienne (élément isoparamétrique)
        J = np.matmul(Bke,coor);
        # Calcul de son inverse
        j = np.linalg.inv(J);
        # Calcul de son déterminant 
        #(utile car l'intégration se fait sur élément de référence plus loin)
        DJ = np.linalg.det(J);
        # Calcul numérique de la matrice Bxy
        Bxy = np.matmul(j,Bke);

        # Construction de la matrice B car élasticité 2D (avec les termes de Bxy)
        # RAPPEL: {epsilon}e=[B]*{U}e
        # [B](3*8) fait le lien entre les 8 déplacements et les 3 déformations
        # Initialisation
        B = np.zeros((3,nnel*2), dtype=float_type);
        # Assignation directe des 16 termes de B non nuls
        B[0,0]=Bxy[0,0];
        B[0,2]=Bxy[0,1];
        B[0,4]=Bxy[0,2];
        B[0,6]=Bxy[0,3];
        B[1,1]=Bxy[1,0];
        B[1,3]=Bxy[1,1];
        B[1,5]=Bxy[1,2];
        B[1,7]=Bxy[1,3];
        B[2,0]=Bxy[1,0];
        B[2,1]=Bxy[0,0];
        B[2,2]=Bxy[1,1];
        B[2,3]=Bxy[0,1];
        B[2,4]=Bxy[1,2];
        B[2,5]=Bxy[0,2];
        B[2,6]=Bxy[1,3];
        B[2,7]=Bxy[0,3];

        # Stockage de la matrice B au point d'intégration en cours
        # Note: cette matrice sera nécessaire pour le calcul des contraintes
        #         gradient(:,:,(ie-1)*(npg)+i)=B(:,:);

        # Contruction numérique de la matrice de rigidite elementaire
        # Note: Sommation des (poids) X (fonction évaluée au point de Gauss)
        # multi_dot fait le produit matriciel de plusiseurs matrices (dans l'ordre)
        ke = ke  + DJ*W[i]*np.linalg.multi_dot([B.transpose(),D,B]);
        me = me + DJ*W[i]*np.tensordot(Nnum.transpose(),Nnum,0);  # matrice masse qui va servir pour le lissage
        # Calcul de la matrice [N](2x8 dans le cas des Q4) à partir de Nnum
        # Initialisation
        N=np.zeros((2,nnel*2), dtype=float_type);
        # Assignation directe des 8 termes de N non nuls
        N[0,0]=Nnum[0];
        N[0,2]=Nnum[1];
        N[0,4]=Nnum[2];
        N[0,6]=Nnum[3];
        N[1,1]=Nnum[0];
        N[1,3]=Nnum[1];
        N[1,5]=Nnum[2];
        N[1,7]=Nnum[3];

        # Calcul du vecteur des forces volumiques 
        # Initialisation
        f=np.zeros(2, dtype=float_type);
        # Calculs des termes du vecteur
        f[0]=rho*fx;
        f[1]=rho*fy;
   
        # Contruction numérique du vecteur de sollicitation élémentaire 
        # Terme des forces volumiques (poids, force magnétique, inertie,...)
        fe  = fe  + DJ*W[i]*N.transpose().dot(f);           
        # Fin de la boucle sur les points d'intégration

    # Assemblage du ke en cours dans KG 
    Kloc=Feeldof(nd,nnel,ddln); # Retourne les # de DDL concernés par l'élément
    KG=Assembl(KG,ke,Kloc);     # Assemble ke dans KG connaissant Kloc
    
    # Assemblage du vecteur sollicitations élémentaire fe
    FG= AssemblF(FG,fe,Kloc);   # Assemble fe dans FG connaissant Kloc
    
    Kloc2=Feeldof(nd,nnel,1); # Retourne les # de DDL concernés par l'élément
    MG=Assembl(MG,me,Kloc2);  # Assemble me dans MG connaissant Kloc2 #matrce masse globale
    #Fin de la boucle sur les éléments intérieurs



# Boucle sur les éléments de contour
# Ajustement des paramètres d'intégration
# On utilise 2 points de Gauss en 1D (intègre exactement ordre 3)
npg=2;  # Nombre de points de Gauss sur les éléments barres
(xlocC, WC) = Gauss(npg)

# Boucle sur les barres de contour
for ie in range(0,nelc):
    # (Car on a mis leurs connectivités à la suite des Q4 dans Connec)   
    # Initialiser à zéro le vecteur élémentaire fe 
    fe = np.zeros(Nuic)

    # Extraire les infos propres à cet élément (coord., connect.,propriétés. ...)
    # Boucle sur les noeuds de l'élément en question 
    nd = np.zeros(nnel, dtype=np.int32)
    coor = np.zeros((nnel,2), dtype=float_type) # Removing z coordinate
    for i in range(0,nnelc):
        nd[i] = mesh.elems_contour[ie,i]
        coor[i,0] = mesh.nodes[nd[i],0]  # Coordonnées des noeuds de l'élément
        coor[i,1] = mesh.nodes[nd[i],1]  # Coordonnées des noeuds de l'élément

    # Longueur de l'élément par le théorème de Pythagore
    Le = sqrt((coor[0,0]-coor[1,0])**2+ (coor[0,1]-coor[1,1])**2)
   
    # Calcul du déterminant de la matrice J
    DJ=Le/2.;

    # Boucle sur les points d'intégration
    for i in range(0,npg):
        # Extraction des coordonnées d'intégration
        k = xlocC[i]; # Ksi du point d'intégration sur la barre en cours
        # Évaluation numérique des fonctions N pour le point et l'élément
        Nnum = Interpolc(k); 

        # Calcul de la matrice [N](2x4 dans le cas des barres) avec Nnum
        # Initialisation
        N=np.zeros((2,nnelc*2), dtype=float_type);
        # Assignation directe des 4 termes de N non nuls
        N[0,0]=Nnum[0];
        N[0,2]=Nnum[1];
        N[1,1]=Nnum[0];
        N[1,3]=Nnum[1];

        # Calcul du vecteur des forces surfaciques
        # Initialisation
        fs=np.zeros(2, dtype=float_type);
        # Calculs des termes du vecteur
        fs[0]=fsx;
        fs[1]=fsy;

        # Contruction numérique du vecteur de charge élémentaire 
        # Il s'agit ici des forces surfaciques (charge répartie, etc.)
        fe = fe  + DJ*WC[i]*N.transpose().dot(fs);
    #Fin de la boucle sur les points d'intégration

    # Assemblage dans le système global aux endroits appropriés
    Kloc = Feeldof(nd,nnelc,ddlnc);
    FG = AssemblF(FG,fe,Kloc);
    #Fin de la boucle sur les éléments de contour

KN = KG
FN = FG

# 10. Traitement des C.L. de déplacement imposé (Dirichlet)
# Imposition des déplacements u  (selon x)
for i in range(0,nclx):          # Boucle sur les noeuds avec u imposé
    nn= CLx_0[i]      # Numéro du noeud en question
    # Calcul de la position de ce DDL
    pos=nn*ddln  # u et v des noeuds avant + 1
    Big= max(abs(KG[pos,:]))
    for j in range(0,ndlt):        # Boucle sur tous les DDL du système
        if(j==pos):
            KG[pos,j]= Big ;     # Met un 1 si c'est le bon DDL
            FG[pos]= Big*CLx_1[i];  # Impose u dans le terme de droite
        else:
            KG[pos,j]=0.0;    # Met un 0 si ce n'est pas le bon DDL

    # Imposition des déplacements v  (selon y)   
    for i in range(0,ncly):            # Boucle sur les noeuds avec v imposé
        nn= CLy_0[i];       # Numéro du noeud en question
        # Calcul de la position de ce DDL
        pos=nn*ddln+1;  # u et v des noeuds avant + 2
        Big= max(abs(KG[pos,:]))
        for j in range(0,ndlt):        # Boucle sur tous les DDL du système
            if(j==pos):
                KG[pos,j]= Big      # Met un 1 si c'est le bon DDL
                FG[pos]= Big*CLy_1[i]  # Impose u dans le terme de droite
            else:
                KG[pos,j]=0.0;    # Met un 0 si ce n'est pas le bon DDL

# 11. RÉSOLUTION DU SYSTÈME EN DÉPLACEMENT (u et v en m par tout)
print_mat(KG/1e11)
sol = np.linalg.solve(KG,FG)

u = np.zeros(nnt, dtype=float_type)
v = np.zeros(nnt, dtype=float_type)
for z in range(0, nnt):
    u[z] = sol[2*z]
    v[z] = sol[2*z+1]

RG = KN.dot(sol)-FN;

# 12. CALCUL DES CONTRAINTES Aus point de Gauss puis on fait un lissage pour les obtenir aux noeuds

GX = np.zeros(nnt)
GY = np.zeros(nnt)
GXY = np.zeros(nnt)

# Ajustement du npg pour les éléments internes
npg=4;

(xloc,W) = Gauss(npg);

# Initialisation de la matrice contenant les contraintes aux nnt noeuds
contraintes_n = np.zeros((3,nnt), dtype=float_type);

for ie in range(0,nelt):  # Pour tous les éléments internes
    # Extraires les déplacements propres à cet élément à partir de ¨sol¨
    # Initialisation du vecteur colonne des déplacements pour cet élément
    uv = np.zeros(Nui)
    gex = np.zeros(nnel)
    gey = np.zeros(nnel)
    gexy = np.zeros(nnel);
    deform = np.zeros(3);
    contraintes = np.zeros(3);
 
    # Repérage de la position des noeuds de cet élément
    nd = mesh.elems[ie,:]  # Numéros des noeuds
    Kloc=Feeldof(nd,nnel,ddln) # Retourne les # de DDL concernés par l'élément

    # extraire le vecteur solution elementaire  des u et v pour les ¨nnel¨ noeuds de l'élément
    for k in range(0,Nui):
        temp=Kloc[k]
        uv[k]=sol[temp]

    nd = np.zeros(nnel, dtype=np.int32)
    coor = np.zeros((nnel,2), dtype=float_type) # Removing z coordinate
    for i in range(0,nnel):
        nd[i] = mesh.elems[ie,i];
        coor[i,0] = mesh.nodes[nd[i],0];  # Coordonnées des noeuds de l'élément
        coor[i,1] = mesh.nodes[nd[i],1];  # Coordonnées des noeuds de l'élément

  
    # ------ Calcul des contraintes aux points d'intégration ----------
    # Boucle sur les points d'intégration
    for ig in range(0,npg):
       
        # Extraction des coordonnees des points d'intégration
        k = xloc[ig,0];  # Ksi du point d'intégration sur le carré en cours
        n = xloc[ig,1];  # Eta du point d'intégration sur le carré en cours
   
        # Évaluation numérique des fonctions N pour le point et l'élément
        Nnum = Interpol(k,n); 
        Bke= BKsiEta(k,n);
   
            
        # Calcul numérique de la matrice Jacobienne (élément isoparamétrique)
        J= np.matmul(Bke,coor)
        # Calcul de son inverse
        j= np.linalg.inv(J)
        # Calcul de son déterminant 
        #(utile car l'intégration se fait sur élément de référence plus loin)
        DJ= np.linalg.det(J);
   
        Bxy = np.matmul(j,Bke);
     
        # Construction de la matrice B car élasticité 2D (avec les termes de Bxy)
        # RAPPEL: {epsilon}e=[B]*{U}e
        # [B](3*8) fait le lien entre les 8 déplacements et les 3 déformations
        # Initialisation
        B=np.zeros((3,nnel*2), dtype=float_type);
        # Assignation directe des 16 termes de B non nuls
        B[0,0]=Bxy[0,0];
        B[0,2]=Bxy[0,1];
        B[0,4]=Bxy[0,2];
        B[0,6]=Bxy[0,3];
        B[1,1]=Bxy[1,0];
        B[1,3]=Bxy[1,1];
        B[1,5]=Bxy[1,2];
        B[1,7]=Bxy[1,3];
        B[2,0]=Bxy[1,0];
        B[2,1]=Bxy[0,0];
        B[2,2]=Bxy[1,1];
        B[2,3]=Bxy[0,1];
        B[2,4]=Bxy[1,2];
        B[2,5]=Bxy[0,2];
        B[2,6]=Bxy[1,3];
        B[2,7]=Bxy[0,3];

        # Calcul et stockage des trois déformations pour le point en cours
        # RAPPEL: {epsilon}e=[B]*{U}e
        deform=B.dot(uv);
        # Calcul et stockage des contraintes au point d'intégration en cours
        # RAPPEL: contraintes=[D]*{epsilon}e
        contraintes=D.dot(deform);  # ICI ON a les contraintes de l'element ie et au point de Gauss IG
        # print(max(contraintes))
        
        gex= gex+W[ig]*Nnum.transpose()*contraintes[0]*DJ;
        gey= gey+ W[ig]*Nnum.transpose()*contraintes[1]*DJ;
        gexy= gexy+W[ig]*Nnum.transpose()*contraintes[2]*DJ;
        
        # fin sur les points d integration
   
    Kloc2=Feeldof(nd,nnel,1); # Retourne les # de DDL concernés par l'élément
     
    GX= AssemblF(GX,gex,Kloc2);   # Assemble gex dans GX connaissant Kloc2
    GY= AssemblF(GY,gey,Kloc2);   # Assemble gey dans GY connaissant Kloc2
    GXY= AssemblF(GXY,gexy,Kloc2);   # Assemble gexy  dans GXY connaissant Kloc2
    # Fin de la boucle sur les éléments
       
# Obtenir les contraintes lissees aux noeuds
sigmax = np.linalg.solve(MG,GX);
sigmay= np.linalg.solve(MG,GY);
sigmaxy= np.linalg.solve(MG,GXY);
contraintes_n[0,:] = sigmax.transpose();
contraintes_n[1,:] = sigmay.transpose();
contraintes_n[2,:] = sigmaxy.transpose();
print("Contrainte maximale sigmax " + str(max(abs(sigmax/10**6))))
print("Contrainte maximale sigmay " + str(max(abs(sigmay/10**6))))
print("Contrainte maximale sigmaxy " + str(max(abs(sigmaxy/10**6))))



# 14. AFFICHAGE DES RÉSULTATS EN DÉPLACEMENT ET CONTRAINTE 
#========================================================= 
export_sol_vtk(mesh,u,v,sigmax/10**6,sigmay/10**6,sigmaxy/10**6)
