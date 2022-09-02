# 4.38 Analytique+schéma implicite+terme source 1D
# Ce programme calcule la solution de l'équation de la chaleur avec un terme source en 1 dimension
# Le graphique renvoyer est l'évolution de la température en fonction de la profondeur au cours du temps
# Il superpose sur un même graphique la solution analytique et la solution à l'aide d'une méthode de différence fini 
# implicite  


using Plots
using SparseArrays
using SpecialFunctions

# Ensemble des fonctions nécessaire au calcule des solutions analytique et avec schéma implicite 

function calcλ1(L,c,DT)
    # Fonction permmettant de determiner la valeur de λ1 dans le cadre du problème 
    # Paramètres d'entrées :
    #   -L: chaleur latente en J kg^-1
    #   -c: capacité thermique en J kg^-1 K^-1
    #   -DT: Différence de température entre le milieu et la surface en Kelvin
    # Paramètre de sortie : 
    #   -λ1 : permet de calculer la valeur de ym [1]  
    λ1_vec=0.1:0.001:2
    fλ1_vec=(exp.(-λ1_vec.^2))./(λ1_vec.*erf.(λ1_vec))

    fλ1=L*sqrt(pi)/(c*DT)
    ind=argmin(abs.(fλ1*ones(length(λ1_vec)).-fλ1_vec))
    λ1=λ1_vec[ind]
    return λ1
end

# Fonction permettant de calculer ym la profondeur maximum de cristalisation pour un temps t prédéfinie 

function calcym(ny,λ1,κ,t_vec,t,nstep)
    # Cette fonction permet de calculer la profondeur ym de cristallion pour un temps t Données
    # Paramètres d'entrées
    #   -ny : nombre de subdivisions de l'intervalle en espace
    #   -λ1 : coefficiant calculer dans la fonction calcλ1
    #   -κ  : coefficient correspondant à k/(ρ*c) [m^2 s^-1]
    #   -t_vec : vecteur temporelle en seconde [s]
    #   t : temps auquel on cherche à calculer ym en années
    #   nstep : nombre de subdivisions de l'intervalle en temps [1]
    # Paramètre de sortie
    #   -ym : profondeur de la solidification en [m]  
    ym=zeros(ny)
    ym_vec=2*λ1*sqrt.(κ*t_vec)
    t_vec=t_vec./(365.25*24*3600)
    time_vec=t*ones(nstep)
    ind=argmin(abs.(time_vec-t_vec))
    ym=ym_vec[ind]
    return ym
end
 
function calcTVA(time,ym,κ,y_vec,λ1,TVA)
    #Cette fonction calcule la solution analytique 
    # Paramètre d'entrée :
    #   -time : temps auquel on cherche à calcule ym en années
    #   -ym: Profondeur de cristallisation en [m]
    #   -κ: Constante thermodynamique (correspond à k/(ρ*Cp))
    #   y_vec: Vecteur des profondeurs [m]
    #   λ1: constante calculer à l'aide de l'abaque
    #   TVA: vecteur des temperatures à l'instant initiale [1]
    #Paramètre de sortie :
    #   -TVA: Vecteur des température [1] selon la profodeur pour un temps t  

    # Calcule de η et de l'indice max de diffusion avec terme source
    η=y_vec./(2*sqrt(κ*time*3600*24*365.25))
    res=argmin(abs.(y_vec.-ym))
    TVA[1:res-1]=erf.(η[1:res-1])./erf(λ1)
    return TVA
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

# Fonction permettant de créer la matrice sparse qui va être utilisé afin de calculer l'équation de diffusion avec
# schéma implicite
function sparse(ny,ρ,k,Cp,dy,dt)
    # On créer la matrice sparse utilisé pour le schéma implicite 
    # Paramètres d'entrée :
    # ny : nombre de subdivision de l'intervalle d'espace
    # ρ : masse volumique [kg m^-3]
    # k : conductivité thermique [kJ m^-1]
    # Cp : capacité thermique massique [kJ kg^-1 K^-1]
    # dy : pas en espace
    # dt :pas de temps 
    # Paramètres de sortie :
    # M : matrice creuse pour résoude le schéma implicite
    μ=zeros(ny-2)
    μ.=k*dt./(Cp.*ρ*dy^2)

    C=[1;ones(ny-2).*(1 .+2*μ);1]
    W=[-ones(ny-2).*μ;0]
    E=[0;-ones(ny-2).*μ]
    M=spdiagm(-1=>W,0=>C,1=>E)
    return M
end

function main()

    # Données:

    # Temperature:
    Tm=1500 #[K] Temperature de la lave
    T0=500 #[K] Temperature à la surface
    DT=Tm-T0 #[K] Différence de température entre le milieu et la surface
    Tr=1500 #[K] Temperature de changement d'état 
    dTr=1 #[K] écart type 

    # Données thermodynamiques:
    L=320e3 #[J kg^-1] 
    c=1e3 #[J kg^-1 K^-1]
    κ=1e-6   #[m^2 s^-1]
    ρ=3000  # masse volumique [kg m^-3] 
    k=ρ*c*κ 

    # Données spatiale 
    ny=401 #Nombre de subdivision de l'intervalle d'espace 
    ymin=0 #[m] Surface 
    ymax=50 #[m] Profondeur max 
    y_vec=LinRange(ymin,ymax,ny) #Vecteur spatiale
    dy=(ymax-ymin)/(ny-1)

    # Données temporelles:
    nstep=201 # Nombre de subdivision de l'intervalle de temps
    time=4
    tmax=time*365.25*24*3600 # Temps finale [s]  
    tmin=0 #Temps initiale [s]
    t_vec=LinRange(tmin,tmax,nstep) # vecteur temporelle 
    dt=(tmax-tmin)/(nstep-1)
    timeb=LinRange(0.1,time,nstep)

    # Condition initiale 
    T_vec=Tm*ones(ny)
    T_vec[1]=T0
    Cp=c*ones(ny)
    dXdT=zeros(ny)
    
    #Calcul de λ1 :
    λ1=calcλ1(L,c,DT)

    for i=1:nstep
        # Calcule de ym et l'indice de y_vec pour lequel y<ym 
        t=timeb[i]
        ym=calcym(ny,λ1,κ,t_vec,t,nstep)
        res=argmin(abs.(y_vec.-ym))
        
        # Calcule de l'avancé de la réaction :
        dXdT=dxdT(T_vec,Tr,dTr)
        
        # Calcule de Cpeff :
        Cp.=c .+L.*dXdT

        #Création de la sparse matrice
        M=sparse(ny,ρ,k,Cp[2:ny-1],dy,dt)
        dropzeros(M)
        
        #Calcule de la solution du schéma implicite et graphe de celle ci 
        Told=copy(T_vec)

        T_vec= M \Told
        
        p2 = plot(T_vec,y_vec,label="implicite", 
                yflip=true,
                title = "Diffusion 1D+terme source", 
                xlabel = "Temperature [K]",
                ylabel = "depth [m]")
        display( plot( p2 ) )
        
        #Calcule de la solution analytique et graphe de celle ci :
        TVA=ones(ny) # Le vecteur est adimensionnée dans le cas de cette méthode (θ=(T-T0)/(Tm-T0))
        TVA[1]=0 # Condition initiale à y=0 T=T0
        TVA=calcTVA(t,ym,κ,y_vec,λ1,TVA)
        TVA=TVA.*DT.+T0 # [K] On passe d'un vecteur de temperature adimensionner à un vecteur  de temperature en Kelvin 
        p1=plot!(TVA,y_vec,label="analytique",
        yflip=true)
        display(plot(p1))


    end
end
main()