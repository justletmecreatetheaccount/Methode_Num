# -------------------------------------------------------------------------
#
# PYTHON for DUMMIES 22-23
# Problème 6
#
# Script de test
#  Vincent Legat
#
# -------------------------------------------------------------------------
#

from numpy import *
import warnings

# ============================================================
# FONCTION A MODIFIER [begin]
#
# u'(t) = 10 * v(t) - 10u(t)
# v'(t) = 28u(t) - v(t) - u(t) * w(t)
# w'(t) = u(t)v(t) - 8/3w(t)
#
#
# Ui+1 = Ui + h/6 * (K1 + 2 * K2 + 2 * K3 + K4)
#
# K1 = f(Ti, Ui)
# K2 = f(Ti + h/2, Ui + h/2 * K1)
# K3 = f(Ti + h/2, Ui + h/2 * K2)
# K4 = f(Ti + h, Ui + h * K3)
#
#

def fn(U):
  F = zeros_like(U)
  F[0] = 10 * U[1] - 10 * U[0]
  F[1] = 28 * U[0] - U[1] - U[0] * U[2]
  F[2] = U[0] * U[1] - (8/3) * U[2]
  return F
f = lambda x: array([10*x[1]-10*x[0],x[0]*(28-x[2])-x[1],x[0]*x[1]-8/3 * x[2]])
def lorenz(Xstart,Xend,Ustart,n):
  X = linspace(Xstart, Xend, n+1)
  h = (Xend - Xstart) / n
  U = zeros((n+1,3), dtype=float64)
  U[0] = asarray(Ustart)
  K1 = zeros(3, dtype=float64)
  K2 = zeros_like(K1)
  K3 = zeros_like(K1)
  K4 = zeros_like(K1)
  #warnings.simplefilter('error')
  U2 = zeros((n+1,3))
  X2 = linspace(Xstart,Xend,n+1)
  U2[0] = Ustart
  for i in range(n):
    K1 = fn (U[i])
    K2 = fn (U[i] + h/2 * K1)
    K3 = fn (U[i] + h/2 * K2)
    K4 = fn (U[i] + h * K3)
    U[i+1] = U[i] + h * (K1 + 2 * K2 + 2 * K3 + K4)/6
    K_1 = f(U2[i])
    K_2 = f(U2[i]+((h/2)*K_1))
    K_3 = f(U2[i]+((h/2)*K_2))
    K_4 = f(U2[i]+(h*K_3))
    U2[i+1] = U[i] + (h * (K_1 + 2*K_2 + 3*K_3 + K_4) /6)

  return X,U,X2,U2


    
#
# FONCTION A MODIFIER [end]
# ============================================================



def main():
  
# ------------------------------------------------------------------------------------ 
#
# Script de test
#
#
# ------------------------------------------------------------------------------------



  from matplotlib import pyplot as plt
  plt.rcParams['axes.facecolor'] = 'lavender'


  plt.figure("Lorenz Equations")
  Xstart = 0; Xend = 100;
  Ustart = [0,1,0]
  n = 10000;

  X,U,X2,U2 = lorenz(Xstart,Xend,Ustart,n)
  plt.plot(U[:,0],U[:,2],'--g',linewidth=0.5)
  plt.plot(U2[:,0],U2[:,2],'.r',linewidth=0.5)
  plt.show()



main()  
