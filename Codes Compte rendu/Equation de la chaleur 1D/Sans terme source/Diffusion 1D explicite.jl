# Diffusion 1D explicite

using Plots

function main()

    # Données:

    # Temperature:
    Tm=1050 #[K] Temperature de la lave
    T0=0 #[K] Temperature à la surface
    DT=Tm-T0 #[K] Différence de température entre le milieu et la surface

    # Données thermodynamiques:
    L=400e3 #[J kg^-1] 
    c=1e3 #[J kg^-1 K^-1]
    κ=0.7*1e-6  #[mm^2 s^-1]
    
    # Données spatiale 
    ny=201 #Nombre de subdivision de l'intervalle d'espace 
    ymin=0 #[m] Surface 
    ymax=30 #[m] Profondeur max 
    y_vec=LinRange(ymin,ymax,ny) #Vecteur spatiale
    dy=(ymax-ymin)/(ny-1)

    # Données temporelles:
    nstep=201 # Nombre de subdivision de l'intervalle de temps
    tmax=4*365.25*24*3600 # Temps finale [s]  
    tmin=0.1*4*365.25*24*3600 #Temps initiale [s]
    t_vec=LinRange(tmin,tmax,nstep) # vecteur temporelle 
    dt=0.249*dy*dy/κ
    
    time =0

    # Condition initiale 
    T_vec=Tm*ones(ny)
    T_vec[1]=T0

    # Boucle temporelle
    T_init=copy(T_vec)
    plot(T_init,y_vec,
    yflip=true)
    for t=2:nstep
        
        time = time + dt
        T_vec[2:ny-1].=((dt*κ)/dy^2).*(T_vec[3:ny].-2*T_vec[2:ny-1].+T_vec[1:ny-2]).+T_vec[2:ny-1]
        p2=plot(T_vec,y_vec,label="T=f(y,t)",
        title="Diffusion 1D explicite; time = $(time/(365.25*24*3600))",
        yflip=true,
        xlabel="Temperature[k]",
        ylabel="depth[m]")
        display(plot(p2))


    end
end
main()