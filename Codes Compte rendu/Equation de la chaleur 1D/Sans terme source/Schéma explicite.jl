#Exercice 4

using Plots

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

x_vec=LinRange(xmin,xmax,nx)
T_vec=Tbg*ones(nx)

for i=1:nx
    if x_vec[i]>=-50 &&  x_vec[i]<=50
        T_vec[i]=Tmax
    end
end
T_init=copy(T_vec)
p=plot(x_vec,T_vec)

for t=1:nstep
    T_vec[2:nx-1].=((dt*kappa)/dx^2).*(T_vec[3:nx].-2*T_vec[2:nx-1].+T_vec[1:nx-2]).+T_vec[2:nx-1] 
    p2=plot(x_vec,T_init)
    p2=plot!(x_vec,T_vec)
    display(plot(p2))
end
