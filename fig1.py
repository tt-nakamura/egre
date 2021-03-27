from sympy import factor
from sympy.abc import x,y
from sympy.vector import CoordSys3D
from sympy.plotting import plot3d
from GaussKH import GaussKH

z = x**3 - 3*x*y**2
e = CoordSys3D('e')
p = x*e.i + y*e.j + z*e.k
K,H = GaussKH(p,x,y)
print('Gauss curvature:', factor(K))
print('mean curvature:', factor(H))

plot3d(z, (x,-1,1), (y,-1,1),
       xlabel='x', ylabel='y', title=r'$z = x^3-3xy^2$')

plot3d(K, (x,-1,1), (y,-1,1),
       xlabel='x', ylabel='y', title='Gauss curvature')
