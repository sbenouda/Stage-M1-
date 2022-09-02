# Ce code renvoie la solution 

using Plots

function main()
    # Donnees

    sigma=100
    sc=sigma*sigma
    Tm=300
    
    nx=101
    ny=51
    nstep=201
    
    kappa=1e-6
    
    xmin=-250
    xmax=250
    ymin=-150
    ymax=150

    dx=(xmax-xmin)/(nx-1)
    dy=(ymax-ymin)/(ny-1)
    dt=0.249*min(dx*dx,dy*dy)/kappa
    # Initialisation vecteurs et matrices
    
    Temp=zeros(ny,nx)

    x_vec=LinRange(xmin,xmax,nx)
    y_vec=LinRange(ymin,ymax,ny)
    xc=x_vec.*x_vec
    yc=y_vec.*y_vec
    for j=1:nx
        Temp[:,j]=Tm*exp.(-(xc[j].+yc)/sc)
    end
    # Temperature initiale

    p1=heatmap(x_vec,y_vec,Temp,c=:jet1,α=0.75,
    title="Temperature a l'instant initiale",
    yflip=false,
    xlabel="Length [Km]",
    ylabel="Width [Km]")
    p1=contour!(x_vec,y_vec,Temp,linestyle=:dash,c=:black,levels=500:500)
    p1=contour!(x_vec,y_vec,Temp,linesstyle=:solid,c=:black,levels=1300:1300)
    display(plot(p1))
    
    
    # Solution analytique

    for t=1:nstep
            for j=1:nx
                Temp[:,j]=Tm/(1+4*t*dt*kappa/sc)*exp.(-(xc[j].+yc)/(sc+4*t*dt*kappa))
            end
        p1=heatmap(x_vec,y_vec,Temp,c=:jet1,α=0.75,
        title="Temperature a l'instant initiale",
        yflip=false,
        xlabel="Length [Km]",
        ylabel="Width [Km]")
        p1=contour!(x_vec,y_vec,Temp,linestyle=:dash,c=:black,levels=500:500)
        p1=contour!(x_vec,y_vec,Temp,linesstyle=:solid,c=:black,levels=1300:1300)
        display(plot(p1))
    end
end

main()