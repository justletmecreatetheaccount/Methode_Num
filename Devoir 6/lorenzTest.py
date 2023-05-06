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

def lorenz(Xstart,Xend,Ustart,n):
  X = linspace(Xstart, Xend, n+1)
  h = (Xend - Xstart) / n
  U = zeros((n+1,3), dtype=float64)
  U[0] = Ustart
  K1 = zeros(3, dtype=float64)
  K2 = zeros_like(K1)
  K3 = zeros_like(K1)
  K4 = zeros_like(K1)
  #warnings.simplefilter('error')
  for i in range(n):
    K1 = fn (U[i])
    K2 = fn (U[i] + h/2 * K1)
    K3 = fn (U[i] + h/2 * K2)
    K4 = fn (U[i] + h * K3)
    U[i+1] = U[i] + (h/6) * (K1 + 2 * K2 + 2 * K3 + K4)
  print(U)
  return X,U


    
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
  from mpl_toolkits.mplot3d.art3d import Line3DCollection
  plt.rcParams['axes.facecolor'] = 'lavender'


  plt.figure("Lorenz Equations")
  Xstart = 0; Xend = 100;
  Ustart = [0,1,0]
  n = 10000;

  X,U = lorenz(Xstart,Xend,Ustart,n)
  plt.plot(U[:,0],U[:,2],'--g',linewidth=0.5)
  plt.show()

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  Xstart = 0; Xend = 100;
  Ustart = [0,1,0]
  n = 10000;
  X, U = lorenz(Xstart,Xend,Ustart,n)
  # Pour avoir un couleur qui s'intensifie de y_0 à y_n
  points = U.reshape(-1,1,3)
  segs = concatenate([points[:-1],points[1:]],axis=1)
  colors = linspace(0, 1, n+1)
  lc = Line3DCollection(segs, cmap=plt.get_cmap('Oranges'))
  lc.set_array(colors)
  ax.add_collection3d(lc)
  x = U[:,0]; y = U[:,1]; z = U[:,2]
  ax.set_xlim(min(x), max(x))
  ax.set_ylim(min(y), max(y))
  ax.set_zlim(min(z), max(z))

  ax.set_xlabel('u(t)')
  ax.set_ylabel('v(t)')
  ax.set_zlabel('w(t)')
  plt.show()




main()  
