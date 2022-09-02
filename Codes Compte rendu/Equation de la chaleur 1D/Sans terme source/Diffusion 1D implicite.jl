# Diffusion 1D implicite 

using Plots
using SparseArrays
function main()

    # Données:

    # Temperature:
    Tm=1050 #[K] Temperature de la lave
    T0=0 #[K] Temperature à la surface
    DT=Tm-T0 #[K] Différence de température entre le milieu et la surface

    # Données thermodynamiques:
    L=400e3 #[J kg^-1] 
    c=1e3 #[J kg^-1 K^-1]
    κ=0.7*1e-6   #[m^2 s^-1]
    
    # Données spatiale 
    ny=201 #Nombre de subdivision de l'intervalle d'espace 
    ymin=0 #[m] Surface 
    ymax=50 #[m] Profondeur max 
    y_vec=LinRange(ymin,ymax,ny) #Vecteur spatiale
    dy=(ymax-ymin)/(ny-1)

    # Données temporelles:
    nstep=201 # Nombre de subdivision de l'intervalle de temps
    tmax=2*365.25*24*3600 # Temps finale [s]  
    tmin=0 #Temps initiale [s]
    t_vec=LinRange(tmin,tmax,nstep) # vecteur temporelle 
    dt=(tmax-tmin)/(nstep-1)
    

    # Condition initiale 
    T_vec=Tm*ones(ny)
    T_vec[1]=T0

    #Création de la sparse matrice
    μ=κ*dt/dy^2
    C=[1;ones(ny-2).*(1+2*μ);1]
    W=[-ones(ny-2).*μ;0]
    E=[0;-ones(ny-2).*μ]

    M=spdiagm(-1=>W,0=>C,1=>E)
    dropzeros(M)

    for t=1:nstep

        Told=copy(T_vec)

        T_vec= M \Told

        p2 = plot(T_vec,y_vec,label="T=f(y,t)", 
                yflip=true,
                title = "Diffusion 1D implicite", 
                xlabel = "Temperature [K]",
                ylabel = "depth [m]")
        display( plot( p2 ) )
    end
end
main()