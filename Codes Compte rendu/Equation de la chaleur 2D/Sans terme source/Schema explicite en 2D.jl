#Schema explicite en 2D
# Il y a une modification a faire dans le code ecrit en premier lieu en 2D,les sx et sy sont a inverser 
#  
using Plots

function main()
    # Donnees

    xmin=-250
    xmax=250
    ymin=-150
    ymax=150

    Tbg=300

    kappa=1e-6

    nx=101
    ny=51
    nstep=101

    dx=(xmax-xmin)/(nx-1)
    dy=(ymax-ymin)/(ny-1)
    dt=0.249*min(dx*dx,dy*dy)/kappa

    sigma=100
    sc=sigma*sigma

    # Initialisation vecteurs et matrices
    
    Temp=zeros(nx,ny)
    x_vec=LinRange(xmin,xmax,nx)
    y_vec=LinRange(ymin,ymax,ny)

    xc=x_vec.*x_vec
    yc=y_vec.*y_vec

    for j=1:ny
        Temp[:,j]=Tbg*exp.(-(xc.+yc[j])/sc)
    end
   

    # Temperature initiale

    p1=heatmap(x_vec,y_vec,Temp',c=:jet1,α=0.75,
    title="Temperature a l'instant initiale",
    yflip=false,
    xlabel="Length [Km]",
    ylabel="Width [Km]")
    p1=contour!(x_vec,y_vec,Temp',linestyle=:dash,c=:black,levels=500:500)
    p1=contour!(x_vec,y_vec,Temp',linesstyle=:solid,c=:black,levels=1300:1300)
    display(plot(p1))

    # Temperature en  fonction du temps

    sx=kappa*dt/(dx*dx)
    sy=kappa*dt/(dy*dy)
    s=1-2*sx-2*sy

    for t=1:nstep
        Tnew=copy(Temp)
        p1=heatmap(x_vec,y_vec,Temp',c=:jet1,α=0.75,
        title="Temperature=f(x,y,t)",
        yflip=false,
        xlabel="Length [Km]",
        ylabel="Width [Km]")
        p1=contour!(x_vec,y_vec,Temp',linestyle=:dash,c=:black,levels=500:500)
        p1=contour!(x_vec,y_vec,Temp',linesstyle=:solid,c=:black,levels=1300:1300)
        display(plot(p1))
        for i=2:nx-1
            for j=2:ny-1
                Tnew[i,j]=sx*(Temp[i+1,j]+Temp[i-1,j])+sy*(Temp[i,j+1]+Temp[i,j-1])+s*Temp[i,j]
            end
        end
        Temp=Tnew
    end
end
main()
    