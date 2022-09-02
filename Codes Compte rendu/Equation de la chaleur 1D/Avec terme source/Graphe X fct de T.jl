# Graphe de X en fonction de la température 

using SpecialFunctions
using Plots

function main()
    dTr=100
    Tr=1050 #[K]
    T=LinRange(0,1500,300)

    ErrT=erfc.(-(T.-Tr)./dTr)

    X1=1.0
    X2=0.0
    dX=X1-X2

    X=X1.-(dX/2)*ErrT

    plot(T,X)
end
main()