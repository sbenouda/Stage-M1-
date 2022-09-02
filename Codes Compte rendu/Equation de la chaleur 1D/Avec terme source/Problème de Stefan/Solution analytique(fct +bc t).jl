# Solution analytique avec des fonctions permmetant de calculer les valeurs de ym et de λ1 puis 
# de les appeler dans le programme principale afin de determiner la temperature selon la profondeur.
# Ce code psoosède une boucle afin de calculer la 
using Plots
using SpecialFunctions

# Fonction permettant de calculer λ1 pour des données thermodynamique défini dans notre problème

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
    #   -κ  : coefficient correspondant à k/(ρ*c) [mm^2 s^-1]
    #   -t_vec : vecteur temporelle en seconde [s]
    #   t : temps auquel on cherche à calculer ym en années
    #   nstep : nombre de subdivisions de l'intervalle en temps [1]
    # Paramètre de sortie
    #   -ym : profondeur de la solidification en [m]  
    ym=zeros(ny)
    ym_vec=2*λ1*sqrt.(κ*t_vec)
    ym_vec=ym_vec*1e-3
    t_vec=t_vec./(365.25*24*3600)
    time_vec=t*ones(nstep)
    ind=argmin(abs.(time_vec-t_vec))
    ym=ym_vec[ind]
    return ym
end

function main()
    # Données 
    #Temperature 
    Tm=1050 #[K] Temperature de la lave
    T0=0 #[K] Temperature en surface 
    DT=Tm-T0 #[K] différence de Temperature entre le milieu et la surface 

    # Données thermodynamiques
    L=400e3 #[J kg^-1] 
    c=1e3 #[J kg^-1 K^-1]
    κ=0.7   #[mm^2 s^-1]
    
    # Données temporelles
    time=4 #Temps en années choisis
    nstep=201 # Nombre de subdivision de l'intervalle de temps
    tmax=time*365.25*24*3600 # Temps finale [s]  
    tmin=0 #Temps initiale [s]
    t_vec=LinRange(tmin,tmax,nstep) # vecteur temporelle 
    timeb=LinRange(0.1,time,nstep) # [yr] Vecteur des temps pour lequelle on calcul la température 

    # Données spatiale 
    ny=201 #Nombre de subdivision de l'intervalle d'espace 
    ymin=0 #[m] Surface 
    ymax=30 #[m] Profondeur max 
    y_vec=LinRange(ymin,ymax,ny) #Vecteur spatiale
    
    # Appel des fonctions calculants λ1 et ym

    λ1=calcλ1(L,c,DT)
    
    #Boucle sur le temps 
    for i=1:nstep
        t=timeb[i]
        ym=calcym(ny,λ1,κ,t_vec,t,nstep)

        #Initialisation du vecteur temperature
        T_vec=ones(ny) # Le vecteur est adimensionnée dans le cas de cette méthode (θ=(T-T0)/(Tm-T0))
        T_vec[1]=0 # Condition initiale à y=0 T=T0

        # Calcule de η et de l'indice max de diffusion avec terme source

        η=y_vec./1e-3./(2*sqrt(κ*t*3600*24*365.25))
        res=argmin(abs.(y_vec.-ym))
        T_vec[1:res-1]=erf.(η[1:res-1])./erf(λ1)

        T_vec=T_vec.*DT.+T0 # [K] On passe d'un vecteur de temperature adimensionner à un vecteur  de temperature en Kelvin 

        p1=plot(T_vec,y_vec,label="T=f(y)",
        title="Temperature en fonction du temps",
        yflip=true,
        xlabel="Temperature [K]",
        ylabel="depth[m]")
        display(plot(p1))
    end
end

main()