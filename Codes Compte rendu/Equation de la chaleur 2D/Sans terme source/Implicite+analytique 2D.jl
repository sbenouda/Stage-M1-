# Test 3 meme graphe implicite et analytique 

using Plots
using SparseArrays

function main()
    # Set up
    # Ici le temps ne depend pas de x et y car schema implicite donc on a 
    # un schema inconditionnelement stable

    xmin=-250
    xmax=250
    ymin=-150
    ymax=150
    sigma=100
    sc=sigma*sigma 
    kappa=1e-6
    Tm=300 #Temperature maximum [C]
    nstep=51
    Ta=365.25*24*3600 #Temps en seconde pour un an


    #Premier graphe pour nx=ny=51


    #Variables que l'on va modifier
    nx=51
    ny=51   
    dx=(xmax-xmin)/(nx-1)
    dy=(ymax-ymin)/(ny-1)
    dt=0.249*min(dx*dx,dy*dy)/kappa
    k=nx*ny

    
    # Initialisation des vecteurs
    x_vec=LinRange(xmin,xmax,nx)
    y_vec=LinRange(ymin,ymax,ny)
    t_vec51=zeros(nstep+1)
    for t=1:nstep+1
        t_vec51[t]=(t-1)*dt
    end

    
    a=ny÷2+1
    b=nx÷2+1

    # Ces deux lignes en commentaires ne verifient que si on est bien au centre. Ne fonctionne que pour des 
    # subdivisions impairs
    # x_vec[a]
    # y_vec[b]

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
        Temp[:,j]=Tm*exp.(-(xc[j].+yc)/sc)
    end
    Tempv=Temp[:,1]
    for j=2:nx
        append!(Tempv,Temp[:,j])
    end


    # Solution analytique

    Tempa=deepcopy(Temp)

    # Vecteurs Temperature au centre pour analytique et schema implicite 

    Tac=zeros(nstep+1) # Temperature au centre solution analytique en fonction du temps  
    Tsic=zeros(nstep+1) # Temperature au centre solution schema implicite en fonction du temps
    
    
    Tac[1]=Tm
    Tsic[1]=Tm


    # Resolution 

    for t=2:nstep+1
        d=1
        Toldv=copy(Tempv)

        Tempv=A\Toldv
         
        # Passage du vecteur Temperature en matrice pour plot

        for j=1:nx
            for i=1:ny
                Temp[i,j]=Tempv[d]
                d=d+1
            end
        end
        for j=1:nx
            Tempa[:,j]=Tm/(1+4*(t-1)*dt*kappa/sc)*exp.(-(xc[j].+yc)/(sc+4*(t-1)*dt*kappa))
        end
        Tac[t]=Tempa[a,b]
        Tsic[t]=Temp[a,b]
    end
    
    plot(t_vec51./Ta,Tac,label="Ta(0,0,t)=f(t)",
        legend=:bottomleft,
        title="Temperature au centre ",
        xlabel="Time in year",
        ylabel="Temperature"
    )
    plot!(t_vec51./Ta,Tsic,linecolor=:red,label="Ts(0,0,t)=f(t),nx=ny=51",markershape=:x)

    # Second graphe pour nx=ny=101


    #Variables que l'on va modifier
    nx=101
    ny=101  
    dx=(xmax-xmin)/(nx-1)
    dy=(ymax-ymin)/(ny-1)
    dt=0.249*min(dx*dx,dy*dy)/kappa
    k=nx*ny
    nstep101=4*nstep
    
    # Initialisation des vecteurs
    x_vec=LinRange(xmin,xmax,nx)
    y_vec=LinRange(ymin,ymax,ny)
    t_vec101=zeros(nstep101+1)
    for t=1:nstep101+1
        t_vec101[t]=(t-1)*dt
    end

    
    a=ny÷2+1
    b=nx÷2+1

    # Ces deux lignes en commentaires ne verifient que si on est bien au centre. Ne fonctionne que pour des 
    # subdivisions impairs
    # x_vec[a]
    # y_vec[b]

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
        Temp[:,j]=Tm*exp.(-(xc[j].+yc)/sc)
    end
    Tempv=Temp[:,1]
    for j=2:nx
        append!(Tempv,Temp[:,j])
    end


    # Solution analytique

    Tempa=deepcopy(Temp)

    # Vecteurs Temperature au centre pour analytique et schema implicite 

    Tacs=zeros(nstep101+1) # Temperature au centre solution analytique en fonction du temps  
    Tsics=zeros(nstep101+1) # Temperature au centre solution schema implicite en fonction du temps
    
    
    Tacs[1]=Tm
    Tsics[1]=Tm


    # Resolution 

    for t=2:nstep101+1
        d=1
        Toldv=copy(Tempv)

        Tempv=A\Toldv
         
        # Passage du vecteur Temperature en matrice pour plot

        for j=1:nx
            for i=1:ny
                Temp[i,j]=Tempv[d]
                d=d+1
            end
        end
        for j=1:nx
            Tempa[:,j]=Tm/(1+4*(t-1)*dt*kappa/sc)*exp.(-(xc[j].+yc)/(sc+4*(t-1)*dt*kappa))
        end
        Tacs[t]=Tempa[a,b]
        Tsics[t]=Temp[a,b]
    end


    plot!(t_vec101./Ta,Tsics,label="Ts(0,0,t)=f(t),nx=ny=101",markershape=:+)

    # Troiseme graphe pour nx=ny=201


    #Variables que l'on va modifier

    nx=201
    ny=201  
    dx=(xmax-xmin)/(nx-1)
    dy=(ymax-ymin)/(ny-1)
    dt=0.249*min(dx*dx,dy*dy)/kappa
    k=nx*ny
    nstep201=16*nstep
    
    # Initialisation des vecteurs
    x_vec=LinRange(xmin,xmax,nx)
    y_vec=LinRange(ymin,ymax,ny)
    t_vec201=zeros(nstep201+1)
    for t=1:nstep201+1
        t_vec201[t]=(t-1)*dt
    end

    
    a=ny÷2+1
    b=nx÷2+1

    # Ces deux lignes en commentaires ne verifient que si on est bien au centre. Ne fonctionne que pour des 
    # subdivisions impairs
    # x_vec[a]
    # y_vec[b]

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
        Temp[:,j]=Tm*exp.(-(xc[j].+yc)/sc)
    end
    Tempv=Temp[:,1]
    for j=2:nx
        append!(Tempv,Temp[:,j])
    end


    # Solution analytique

    Tempa=deepcopy(Temp)

    # Vecteurs Temperature au centre pour analytique et schema implicite 

    Tact=zeros(nstep201+1) # Temperature au centre solution analytique en fonction du temps  
    Tsict=zeros(nstep201+1) # Temperature au centre solution schema implicite en fonction du temps
    
    
    Tact[1]=Tm
    Tsict[1]=Tm


    # Resolution 

    for t=2:nstep201+1
        d=1
        Toldv=copy(Tempv)

        Tempv=A\Toldv
         
        # Passage du vecteur Temperature en matrice pour plot

        for j=1:nx
            for i=1:ny
                Temp[i,j]=Tempv[d]
                d=d+1
            end
        end
        for j=1:nx
            Tempa[:,j]=Tm/(1+4*(t-1)*dt*kappa/sc)*exp.(-(xc[j].+yc)/(sc+4*(t-1)*dt*kappa))
        end
        Tact[t]=Tempa[a,b]
        Tsict[t]=Temp[a,b]
    end

    plot!(t_vec201./Ta,Tsict,label="Ts(0,0,t)=f(t),nx=ny=201",markershape=:o)


# Quatrieme graphe pour nx=ny=401


    #Variables que l'on va modifier
    nx=401
    ny=401  
    dx=(xmax-xmin)/(nx-1)
    dy=(ymax-ymin)/(ny-1)
    dt=0.249*min(dx*dx,dy*dy)/kappa
    k=nx*ny
    nstep401=64*nstep
    
    # Initialisation des vecteurs
    x_vec=LinRange(xmin,xmax,nx)
    y_vec=LinRange(ymin,ymax,ny)
    t_vec401=zeros(nstep401+1)
    for t=1:nstep401+1
        t_vec401[t]=(t-1)*dt
    end

    
    a=ny÷2+1
    b=nx÷2+1

    # Ces deux lignes en commentaires ne verifient que si on est bien au centre. Ne fonctionne que pour des 
    # subdivisions impairs
    # x_vec[a]
    # y_vec[b]

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
        Temp[:,j]=Tm*exp.(-(xc[j].+yc)/sc)
    end
    Tempv=Temp[:,1]
    for j=2:nx
        append!(Tempv,Temp[:,j])
    end


    # Solution analytique

    Tempa=deepcopy(Temp)

    # Vecteurs Temperature au centre pour analytique et schema implicite 

    Tacq=zeros(nstep401+1) # Temperature au centre solution analytique en fonction du temps  
    Tsicq=zeros(nstep401+1) # Temperature au centre solution schema implicite en fonction du temps
    
    
    Tacq[1]=Tm
    Tsicq[1]=Tm


    # Resolution 

    for t=2:nstep401+1
        d=1
        Toldv=copy(Tempv)

        Tempv=A\Toldv
         
        # Passage du vecteur Temperature en matrice pour plot

        for j=1:nx
            for i=1:ny
                Temp[i,j]=Tempv[d]
                d=d+1
            end
        end
        for j=1:nx
            Tempa[:,j]=Tm/(1+4*(t-1)*dt*kappa/sc)*exp.(-(xc[j].+yc)/(sc+4*(t-1)*dt*kappa))
        end
        Tacq[t]=Tempa[a,b]
        Tsicq[t]=Temp[a,b]
    end


    plot!(t_vec401./Ta,Tsicq,label="Ts(0,0,t)=f(t),nx=ny=401")


end
main()