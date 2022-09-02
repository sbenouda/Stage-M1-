# Schema implicite 2D+CL unsulating +chaleur latente 

# Ce programme permet de résoudre l'équation de la chaleur sans second membre en 2 dimensions avec CL sans flux 
# Il fonctionne egalement dans le cas oú il n'y a pas de chaleur latente mais il n'est pas conseille d'utiliser cette version dans ce cas 
# car elle est plus longue dû au fait que la matrice creuse utiliser pour resoudre l'equation est appele dans la boucle de temps et 
# cela rallonge considerablememt le temps de calcul selon la dimension des subdivions selon l'axe des x et celui des y
 
# La fonction Cpeff n'est pas utilise ici

# Lorsque L=0 (pas de chaleur latente) le code fonctionne, il ne fonctionne pas 
# lorsque L est different de 0 


using Plots
using SparseArrays

function sparse(nx,ny,k,sx,sy,s)
    # Fonction qui va créer la matrice creuse utiliser afin de calculer la solution de l'équation de la chaleur en 2D
    # Paramètres d'entrées :
        # nx : nombre de subdivisons de l'intervalle en espace selon l'axe des x (horizontale)
        # nx : nombre de subdivion de l'intervalle en espace selon l'axe des y (correpond à la verticale)
        # k : correspond aux produits des deux subdivisions en espace selon l'axe des x et l'axe des y_vec
        # sx : constante devant les termes dependant de l'axe des x dans le schéma implicite 
        # sy : constante devant les termes dependant de l'axe des y dans le schéma implicte 
        # s : constante devant les termes depandant de l'axe des y et l'axe des x dans le schéma implicite 
    # Parametres de sorties :
        # A : Matrice creuse regroupant l'ensemble des constantes devant chaque terme du schéma implicite

    # Creation des vecteurs 

    C=vcat(ones(ny+1),s[ny+2:(nx-1)*ny-1].*ones(k-2*(ny+1)),ones(ny+1))
    N=vcat(zeros(ny),-sy[ny+1:(nx-1)*ny-2].*ones(k-2*(ny+1)),zeros(ny+1))
    S=vcat(zeros(ny+1),-sy[ny+2:(nx-1)*ny-1].*ones(k-2*(ny+1)),zeros(ny))
    W=vcat(0,-sx[2:k-ny].*ones(k-ny-1))
    E=vcat(-sx[1:k-ny-1].*ones(k-ny-1),0)

    # Modification des vecteurs pour avoir la matrice A souhaite
    for l=2:nx-2
        # Termes en ij
        C[l*ny]=1
        C[l*ny+1]=1

        # Termes en i-1,j
        N[l*ny-1]=0
        N[l*ny]=0

        # Termes en i+1,j
        S[l*ny]=0
        S[l*ny+1]=0

        # Termes en i,j-1
        W[l*ny-ny]=0
        W[l*ny+1-ny]=0

        # Termes en i,j+1
        E[(l-1)*ny+ny]=0    
        E[(l-1)*ny+1+ny]=0
    end

    # Modification sur W et E

    W[k-2*ny:k-ny]=zeros(length(W[k-2*ny:k-ny]))   
    E[1:ny+1]=zeros(length(E[1:ny+1]))                  

    # Creation de la marice A
    A=spdiagm(-ny=>W,-1=>N,0=>C,1=>S,ny=>E)
end

function switch(opt,ny,nx,A,sx,sy,s)
    # Fonction qui fait office de switch selon si il y a une CL unsulating ou non 
    # Arguments d'entrées :
        # opt : Entier qui selon la valeur nous dit si on a une CL unsulating sur la verticale, l'horizontale 
        #       ou si il n'y a pas de CL unsulating 
        # ny : Nombre de subdivision de l'intervalle en espace pour la partie selon l'axe des y
        # nx : Nombre de subdivion de l'intervalle en espace pour la partie selon l'axe des x 
        # A : Matrice creuse créer pour résoudre l'équation avec un schéma implicite
        # sx : Constante devant les termes temperatures a cote en horizontale du terme recherche 
        # sy : Constante devant les termes temperatures a cote en verticale du terme recherche
        # s : Constante devant le terme temperature a l'instant suivant dans le temps
    # Arguments de sorties :
        # Modifie les termes sur les bords si CL ou ne fait rien selon le switch de la matrice creuse A  

    # Dans le cas oú on a une CL unsulating sur les bords gauche et droit de la matrice des temperatures
    if opt==1
        # Bord gauche avec CL sans flux 
        for l=2:ny-1
            A[l,l]=s[l]
            A[l,l-1]=-sy[l-1]
            A[l,l+1]=-sy[l+1]
            A[l,l+ny]=-2*sx[l+ny]
        end
        # Bord droit avec Cl sans flux 
        for k=(nx-1)*ny+2:ny*nx-1
            A[k,k]=s[k]
            A[k,k-1]=-sy[k-1]
            A[k,k+1]=-sy[k+1]
            A[k,k-ny]=-2*sx[k-ny]
        end
    end
        # Dans le cas ou on une CL sans flux sur les bords nord et sud  
    if opt==2
        for l=2:nx-1
            # Bord du bas avec une Cl sans flux 
            A[l*ny,l*ny]=s[l*ny] # Terme sur la diagonale centrale, correspond a la temperature a l'endroit ou l'on calcule  
            A[l*ny,l*ny-ny]=-sx[l*ny-ny] # Terme sur l'extra diagonale a gauche, correspond a la temperature sur le point a gauche de la croix 
            A[l*ny,l*ny+ny]=-sx[l*ny+ny] # Terme sur l'extra diagonale a droite, correspond a la temperature sur le point a droite de la croix  
            A[l*ny,l*ny-1]=-2*sy[l*ny-1] # Terme sur la diagonale du dessous de celle centrale, correpons a la temperature sur le point nord de la croix  
            
            # Bord du dessus avec une CL sans flux  
            # On prend l-1 car les termes sur le bord du dessus dans le numbering sont classe de nx+1,2nx+1 a (ny-2)nx+1 (on ne prend pas en compte les coins)
            A[(l-1)*ny+1,(l-1)*ny+1]=s[(l-1)*ny+1] # Terme sur la diagonale centrale, correspond a la temperature a l'endroit ou l'on calcule 
            A[(l-1)*ny+1,(l-1)*ny+1-ny]=-sx[(l-1)*ny+1-ny] # Terme sur l'extra diagonale a gauche, correspond a la temperature sur le point a gauche de la croix 
            A[(l-1)*ny+1,(l-1)*ny+1+ny]=-sx[(l-1)*ny+1+ny] # Terme sur l'extra diagonale a droite, correspond a la temperature sur le point a droite de la croix 
            A[(l-1)*ny+1,(l-1)*ny+1+1]=-2*sy[(l-1)*ny+1+1] # Terme sur la diagonale du dessous de celle centrale, correpons a la temperature sur le point nord de la croix
        end
    end
end 

function dxdT(T,Tr,dTr)
    # Cette fonction calcule la dérivée de l'avancé du 
    # changement d'état au cours du temps 
    # Paramètre d'entrée:
    # T :Température au cours du temps 
    # Paramètre de sortie :
    # dXdT: Avancé de la réaction(compris entre 0 et 1)
    
    # Données dans le cas de la lave 
    dXdT=1.0*exp.(-(-T.+Tr).^2 ./dTr.^2)./(sqrt(pi)*dTr)
    return dXdT
end

function Cpeff(c,L,dxdt,ny)
    # Cette fonction doit calculer la capacite thermique en chaque point de notre modele 
    # Parametres d'entrees :
    #   c : capacite thermique en J kgˆ-1 Kˆ-1
    #   L : chaleur latente en J kgˆ-1
    #   dxdt : matrice de taille ny,nx regroupant l'avance de la reaction en chacun des points 
    #   nx : nombre de subdivision de l'axe des x (largeur)
    #   ny : nombre de subdivion de l'axe des y  (profondeur)
    # Parametre de sortie :
    #   Cp : Vecteur regroupant l'ensemble des capactites thermiques selon la profondeur a un point de l'axe des x 
    
    # Initialisation de la matrice Cp
    Cp=c*ones(ny)
    # Calcule 
    Cp.= c.+L*dxdt
    return Cp
end

function main()
    # Données:

    # Temperature:
    Tm=1500 #[K] Temperature de la lave
    T0=500 #[K] Temperature à la surface
    DT=Tm-T0 #[K] Différence de température entre le milieu et la surface
    Tr=1500
    dTr=1

    # Données thermodynamiques:
    L=320e3 #[J kg^-1] 
    #L=0 # [J kg^-1] Si il n'y pas de terme source
    c=1e3 #[J kg^-1 K^-1]
    κ=1e-6   #[m^2 s^-1]
    ρ=3000
    k=ρ*c*κ

    # Données spatiale
    # depth 
    ny=201 #Nombre de subdivision de l'intervalle d'espace 
    ymin=0 #[m] Surface 
    ymax=50 #[m] Profondeur max 
    y_vec=LinRange(ymin,ymax,ny) #Vecteur spatiale
    dy=(ymax-ymin)/(ny-1)

    # Largeur 
    nx=201
    xmin=0
    xmax=50
    x_vec=LinRange(xmin,xmax,nx)
    dx=(xmax-xmin)/(nx-1)

    kd=nx*ny # Nombre de lignes et de colonnes de la matrice creuse permettant de resoudre l'equation 

    # Centre de la matrice temperature (cas impaire)
    a=ny÷2+1    # a represente la ligne qui caracterise la ligne au centre pour la matrice temperature
    b=nx÷2+1    # b represente la colonne qui caracterise la colonne au centre pour la maatrice temperature    

    # Données temporelles:
    nstep=201 # Nombre de subdivision de l'intervalle de temps
    time=4
    tmax=time*365.25*24*3600 # Temps finale [s]  
    tmin=0 #Temps initiale [s]
    t_vec=LinRange(tmin,tmax,nstep) # vecteur temporelle 
    dt=(tmax-tmin)/(nstep-1)
    timeb=LinRange(0.1,time,nstep)

    # CL sans flux
    # Si opt=0 il n'y a pas de CL sans flux 
    # Si opt=1  cela signifie que l'on choisis des CL sans flux sur les bords est et ouest de notre modele  
    # Si opt=2 cela signifie que l'on choisis des CL sans flux sur les bords nord et sud de notre modele 
    opt=2
    # Initialisation des vecteurs et de la matrice Temperature
    Temp=Tm*ones(ny,nx)
    # Cas ou opt=1 
    if opt==1
        Temp[1,:].=T0
        Temp[ny,:].=Tm
    end

    # Cas ou opt=2
    if opt==2
        Temp[:,1].=T0
        Temp[:,nx].=Tm
    end

    # Vecteur temperatures
    Tempv=Temp[:,1]
    for i=2:nx
        append!(Tempv,Temp[:,i])
    end

    # Initialisations des constantes sx,sy et s contenuent dans le modele 
    sx=zeros(ny,nx)
    sy=zeros(ny,nx)
    s=zeros(ny,nx)

    #Initialisation de la matrice avancee de la reaction 
    dXdT=zeros(ny,nx)

    # Initialisation de la matrice contenant les capacite thermiques
    Cp=c*ones(ny,nx)

    # Partie Resolution
    for t=1:nstep
        d=1 # Compteur pour passer du vecteur temperature a la matrice 

        # Avance de la reaction a chaque instant t et modification de Cp
        for i=1:nx
            dXdT[:,i]=dxdT(Temp[:,i],Tr,dTr)
            Cp[:,i]=c.+L*dXdT[:,i]
            sx[:,i]=k*dt./(ρ.*Cp[:,i]*dx^2)
            sy[:,i]=k*dt./(ρ.*Cp[:,i]*dy^2)
            s[:,i]=1 .+2*sx[:,i].+2*sy[:,i]
        end

        # Passage de matrice a vecteur des elements s,sx et sy 
        # On fait ce passage afin que cela soit plus facile dans les indices pour creer la matrice creuse
        sxv=sx[:,1] # Initialisation du vecteur sx de taille nx*ny
        syv=sy[:,1] # Initialisation du vecteur sy de taille nx*ny
        sv=s[:,1] # Initialisation du vecteur s de taille nx*ny
        for j=2:nx
            sxv=append!(sxv,sx[:,j])
            syv=append!(syv,sy[:,j])
            sv=append!(sv,s[:,j])
        end

        # Creation de la matrice creuse utilise afin de resoudre l'equation de la chaleur en implicite 
        A=sparse(nx,ny,kd,sxv,syv,sv)
        dropzeros(A)
        

        # Prise en compte des CL sans flux 
        switch(opt,ny,nx,A,sxv,syv,sv)
        
        # Resolution
        Toldv=copy(Tempv)

        Tempv=A\Toldv

        # Passage du vecteur Temperature en Matrice Temperature
        for j=1:nx
            for i=1:ny
                Temp[i,j]=Tempv[d]
                d=d+1 
            end
        end
        p1=heatmap(x_vec,y_vec,Temp,c=:jet1,α=0.75,
        title="Equation de la chaleur 2D en implicite",
        yflip=true,
        xlabel="Length[km]",
        ylabel="depth [km]")
        p1=contour!(x_vec,y_vec,Temp,linestyle=:dash,c=:black,levels=500:500)
        display(plot(p1))


        #
        #if opt==1
        #  Temp_vec=Temp[:,b]
        #  p1=plot(Temp_vec,y_vec,label="T=f(x,y,t)",
        #    yflip=true,
        #    title="Diffusion avec terme source et CL sans flux ",
        #    xlabel="Temperature [K]",
        #    ylabel="depth[m]")
        #    display(plot(p1))
        #end

        #if opt==2
        #    Temp_vec=Temp[a,:]
        #   p1=plot(x_vec,Temp_vec,label="T=f(x,y,t)",
        #    title="Diffusion avec terme source et Cl sans flux",
        #    xlabel="largeur [m]",
        #    ylabel="Temperature[K]")
        #    display(plot(p1))
        #end
    end
end
main()
