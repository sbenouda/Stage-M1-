# Figure 4.33

using Plots
using SpecialFunctions
function main()

    DT=1050 #[K] Tm-T0

    L=400e3 #[J kg^-1]
    c=1e3 #[J kg^-1 K^-1]
    κ=0.7   #[mm^2 s^-1]
    nstep=201
    ny=201

    # Abaque pour deteminer λ1

    λ1_vec=0.1:0.001:2
    fλ1_vec=(exp.(-λ1_vec.^2))./(λ1_vec.*erf.(λ1_vec))
    plot(λ1_vec,fλ1_vec,label="f(λ1)",
    title="Abaque λ1",
    xlabel="λ1",
    ylabel="e(-λ1^2))/(λ1*erf(λ1))",
    #xlims=(0.9,1),
    #ylims=(0,5)
    )

    # Determination λ1 dans cette Exercice

    fλ1=L*sqrt(pi)/(c*(DT))
    plot!(λ1_vec,fλ1*ones(length(λ1_vec)))

    ind=argmin(abs.(fλ1*ones(length(λ1_vec)).-fλ1_vec))
    λ1=λ1_vec[ind]

    # Résolution 

    Tmax=4*365.25*24*3600 #temps max 4 ans en s 
    Tmin=0
    t_vec=LinRange(Tmin,Tmax,nstep)

    ymin=0 #[m]
    ymax=50 #[m]
    y_vec=LinRange(ymin,ymax,ny)
    T_vec=ones(ny)
    T_vec[1]=0
    ym_vec=zeros(ny)
    for t=1:nstep
        ym_vec[t]=2*λ1*sqrt(κ*t_vec[t])
    end 

    ym_vec=ym_vec*10^-3 #Passage de mm en m
    t_vec=t_vec./(365.25*24*3600)

    plot(t_vec,ym_vec,label="ym=f(t)",
    title="Thiknesses of the solidifying crust on the lava lake",
    yflip=true,
    xlabel="t(yr)",
    ylabel="ym(m)")
end
main()