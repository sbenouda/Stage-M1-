

using Plots
using SparseArrays

function main()
    #Set up
    xmin=-250
    xmax=250
    nx=101
    Tbg=300 #[°C]
    Tmax=1000 #[°C]
    kappa=1e-6
    nstep=101
    dx=(xmax-xmin)/(nx-1)
    dt=0.249*dx*dx/kappa

    #Initialisation des vecteurs
    x_vec=LinRange(xmin,xmax,nx)
    T_vec=Tbg*ones(nx)

    #Initialisation du vecteur initiale
    for i=1:nx
        if x_vec[i]>=-50 &&  x_vec[i]<=50
            T_vec[i]=Tmax
        end
    end

    # plot des conditions initiale
    p = plot(x_vec, T_vec, label="Tinit = f(x)", 
    marker = 2,
    title = "Temperature", 
    xlabel = "x [m]",
    ylabel = "T [°C]")

    Temp_init = copy(T_vec)
    #Création de la sparse matrice
    mu=kappa*dt/dx^2
    C=[1;ones(nx-2).*(1+2*mu);1]
    W=[-ones(nx-2).*mu;0]
    E=[0;-ones(nx-2).*mu]

    M=spdiagm(-1=>W,0=>C,1=>E)
    dropzeros(M)

    for t=1:nstep

        Told=copy(T_vec)

        T_vec= M \Told

        p2 = plot(x_vec, Temp_init, 
                marker = 2,
                title = "Temperature", 
                xlabel = "x [m]",
                ylabel = "T [°C]")
            p2 = plot!(x_vec, T_vec, marker = 1)
            display( plot( p2 ) )
    end
end 

main()