#
# PYTHON for DUMMIES 22-23
# Problème 3
#
# Splines cubiques périodiques
#
#  Vincent Legat
#
# -------------------------------------------------------------------------
# 


from numpy import *
from numpy.linalg import solve


# ============================================================
# FONCTIONS A MODIFIER [begin]
#
#

def spline(x,h,U):

#
# -0- Construction des abscisses (y-compris celle qui correspond au retour au point de depart)
#     Exemple U = [U_0,U_1,U_2] => n=3 => X = [0,h,2h,3h]
#   
#
  n = size(U)
  X = arange(0,n+1)*h


#
# -1- Detection de l'intervalle dans lequel un élément de x appartient :-)
#     Vous pouvez modifier cette partie de la fonction, il y a des versions plus efficaces :-)
#     Certains étudiants m'ont parlé de trucs rigolos super fast quand il y a beaucoup de points
#     Et Nathan m'a trouvé un truc super rapide et méga cryptique :-)
#  
#     Mais n'oubliez pas : il faut étudier et travailler l'analyse et pas passer trop
#     de temps à des problèmes python qui ne rapporte que 0.2 points au final : think about it !
#

  i = zeros(len(x),dtype=int)
  for j in range(1,n):
      i[X[j]<=x] = j
      
#
# A MODIFIER ..... [begin]
# Ici, on renvoie une interpolation linéaire par morceaux :-)
# Un peu trop simple non !
#

  matrix = eye(n,n)*4 + eye(n,n, k=-1) + eye(n,n, k=1)
  matrix[0][n-1] = 1
  matrix[n-1][0] = 1
  matrix = matrix * h**2/6
  

  vector = []
  for k in range(n-1):
    vector.append(U[k-1] - 2*U[k] + U[k+1])
  
  vector.append(U[n-2] - 2*U[n-1] + U[0])
  ddU = solve(matrix, vector)
  ddU = append(ddU, ddU[0])
  print(vector)
  U = append(U, U[0])

  return (ddU[i] * (X[i + 1] - x) ** 3) / (6 * h) + (ddU[i + 1] * (x - X[i]) ** 3) / (6 * h) + ((U[i]/h) - (ddU[i] * (h / 6))) * (X[i + 1] - x) + ((U[i + 1]/h) - (ddU[i+1] * (h / 6))) * (x - X[i])
         
#
# A MODIFIER ..... [end]
#

#
# FONCTIONS A MODIFIER [end]
# ============================================================
#
# -1- Interpolation d'un cercle :-)     
#

def main() :

  from matplotlib import pyplot as plt
  plt.rcParams['toolbar'] = 'None'
  plt.rcParams['figure.facecolor'] = 'lavender'

  n = 4;
  h = 3*pi/(2*(n+1));
  T = arange(0,3*pi/2,h)
  X = cos(T); Y = sin(T)

  fig = plt.figure("Splines cubiques et cercle :-)")
  plt.plot(X,Y,'.r',markersize=10)
  t = linspace(0,2*pi,100)
  plt.plot(cos(t),sin(t),'--r')

  t = linspace(0,3*pi/2,100)
  plt.plot(spline(t,h,X),spline(t,h,Y),'-b')
  plt.axis("equal"); plt.axis("off")
  plt.show()
 
if __name__ == '__main__': 
  main()
