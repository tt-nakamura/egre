from sympy import sin,cos,cosh,sinh,pi
from sympy.abc import u,v
from sympy.vector import CoordSys3D
from sympy.plotting import plot3d_parametric_surface
from GaussKH import GaussKH

e = CoordSys3D('e')

# catenoid
x,y,z = cosh(u)*cos(v), cosh(u)*sin(v), u
p = x*e.i + y*e.j + z*e.k
K,H = GaussKH(p)
print('Gauss curvature:', K)
print('mean curvature:', H)
plot3d_parametric_surface(
    x,y,z,(u,-2,2),(v,-pi,pi),
    xlabel='x', ylabel='y', title='catenoid')

# helicoid
x,y,z = sinh(u)*cos(v), sinh(u)*sin(v), v
p = x*e.i + y*e.j + z*e.k
K,H = GaussKH(p)
print('Gauss curvature:', K)
print('mean curvature', H)
plot3d_parametric_surface(
    x,y,z,(u,-2,2),(v,-pi,pi),
    xlabel='x', ylabel='y', title='helicoid')
