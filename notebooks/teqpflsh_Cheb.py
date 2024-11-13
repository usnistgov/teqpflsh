import ChebTools
import teqpflsh
import numpy as np 
import matplotlib.pyplot as plt

f = lambda x: np.exp(x)
f = lambda x: x**2-0.5

f = lambda x: np.cos(x)
xmin, xmax = [-10, 10]

CTce = ChebTools.generate_Chebyshev_expansion(120, f, xmin, xmax)

ces = ChebTools.dyadic_splitting(12, f, xmin, xmax, 3, 1e-12, 12)
ces_ = [teqpflsh.ChebyshevExpansion(xmin=_.xmin(), xmax=_.xmax(), coeff=_.coef()) for _ in ces]
for _ in ces_:
    print(_.xmin, _.xmax)
ca = teqpflsh.ChebyshevApproximation(expansions=ces_)
