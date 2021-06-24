# Grace, 2019/02/18

# reference:
#   YouTube lecture: https://www.youtube.com/watch?v=RHkWFlIhNHo&index=3&list=PL1c_abw8OtQ_U3fuvA8UEB8Kq4r2lIsb3
#   coefficient: https://aip.scitation.org/doi/pdf/10.1063/1.1677028?class=pdf
#   Gaussian basis: https://www.worldscientific.com/doi/abs/10.1142/9789812832115_0001
import numpy as np
from constant import *


def Gen_SOrbital(Ax, Ay, Az, alpha):
    # Build_SOrbital, single primitive Gaussian orbital
    N = (2 * alpha / pi) ** 0.75
    return [Ax, Ay, Az, alpha, N]


def Gen_Basis(Z, AL):
    # Build_Basis
    NAtoms = len(AL)
    N = 0
    nb = 0
    for index in range(NAtoms):
        x0 = AL[index][0]
        y0 = AL[index][1]
        z0 = AL[index][2]
        N = N + 1
        S = [[18.7311370, 0.0334946],
             [2.8253937, 0.2347269],
             [0.6401217, 0.81375733]]
        nb = nb + 1
        basis[nb].Gen_SOrbital
    return basis, N


Build_SOrbital = []
Build_SOrbital.extend([Gen_SOrbital(0, 0, -1, 7)])
Z = [1, 1]
AL = [[0, 0, 0], [0, 0, 3]]
Build_Basis = []
Build_Basis.extend([Gen_Basis(Z, AL)])

# Gen_SOrbital(1, 2, 3, 4)
print(Build_SOrbital[0][4])
