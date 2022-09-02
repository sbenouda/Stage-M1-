from sympy import *
import numpy as np
init_printing()

dTr, dA, dn, dQ, Tr, T  = symbols('dTr, dA, dn, dQ, Tr, T')

#########
# X et ses derivées w.r.t pressures - WITHOUT KINETICS


ErrT     = erfc( - (T - Tr) / dTr )
X1       = 0.0;
X2       = 1.0;
dX       = X1-X2;
X        = X1 - dX/2 * ErrT;
dXdT     = X.diff(T)

print( 'dXdP   = ' + octave_code(dXdT) + ';')
print( 'dXdP   = ' + ccode(dXdT) + ';')

