from sympy import symbols,Function
from sympy.abc import u,v
from sympy.vector import CoordSys3D
from GaussKH import GaussKH
from Brioschi import metric,Brioschi

x,y,z = symbols('p,y,z', cls=Function)
e = CoordSys3D('e')
p = x(u,v)*e.i + y(u,v)*e.j + z(u,v)*e.k # any surface
K,H = GaussKH(p)
K1 = Brioschi(metric(p))
print((K-K1).simplify()) # 0 if Theorema Egregium holds
