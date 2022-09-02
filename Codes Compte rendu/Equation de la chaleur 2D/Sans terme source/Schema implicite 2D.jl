# Fonction permettant de resoudre l'equation de la chaleur en 2D  a l'aide d'un schema implicite
# Ceci est un cas particulier, ce programme n'est pas le plus aboutit et ne permet de resoudre l'equation de la chaleur uniquement
# dans le cas ou la capacite thermique c est une constante

using Plots
using SparseArrays

function main()
    # Set up
    xmin=-250
    xmax=250
    ymin=-150
    ymax=150
    nx=101
    ny=51
    Tbg=300
    kappa=1e-6
    nstep=201
    dx=(xmax-xmin)/(nx-1)
    dy=(ymax-ymin)/(ny-1)
    dt=0.249*min(dx*dx,dy*dy)/kappa
    k=nx*ny
    sigma=100
    sc=sigma*sigma

    # Initialisation des vecteurs
    x_vec=LinRange(xmin,xmax,nx)
    y_vec=LinRange(ymin,ymax,ny)

    sx=kappa*dt/(dx*dx)
    sy=kappa*dt/(dy*dy)
    
    s=1+2*sx+2*sy

    # Creation de la sparse matrice
    # Creation des vecteurs 

    C=vcat(ones(ny+1),s*ones(k-2*(ny+1)),ones(ny+1))
    N=vcat(zeros(ny),-sy*ones(k-2*(ny+1)),zeros(ny+1))
    S=vcat(zeros(ny+1),-sy*ones(k-2*(ny+1)),zeros(ny))
    W=vcat(0,-sx*ones(k-ny-1))
    E=vcat(-sx*ones(k-ny-1),0)

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
    dropzeros(A)

    #Condition initiale et affichage de celle ci
    
    Temp=zeros(ny,nx)
    #Tempv=zeros(k,1)
    xc=x_vec.*x_vec
    yc=y_vec.*y_vec

    for j=1:nx
        Temp[:,j]=Tbg*exp.(-(xc[j].+yc)/sc)
    end
    Tempv=Temp[:,1]
    for j=2:nx
        append!(Tempv,Temp[:,j])
    end


    p1=heatmap(x_vec,y_vec,Temp,c=:jet1,α=0.75,
    title="Equation de la chaleur 2D en implicite",
    yflip=false,
    xlabel="Length[km]",
    ylabel="Width [km]")
    p1=contour!(x_vec,y_vec,Temp,linestyle=:dash,c=:black,levels=500:500)
    p1=contour!(x_vec,y_vec,Temp,linestyle=:solid,c=:black,levels=1300:1300)
    display(plot(p1))


    # Resolution

    for t=1:nstep
        d=1
        Toldv=copy(Tempv)

        Tempv=A\Toldv
         
        # Passage du vecteur Te;perature en matrice pour plot

        for j=1:nx
            for i=1:ny
                Temp[i,j]=Tempv[d]
                d=d+1
            end
        end
        p1=heatmap(x_vec,y_vec,Temp,c=:jet1,α=0.75,
        title="Equation de la chaleur 2D en implicite",
        yflip=false,
        xlabel="Length[km]",
        ylabel="depth [km]")
        p1=contour!(x_vec,y_vec,Temp,linestyle=:dash,c=:black,levels=500:500)
        display(plot(p1))
    end
end
    

main()