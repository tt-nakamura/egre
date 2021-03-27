from sympy import diff,simplify
from sympy.abc import u,v

def GaussKH(p,u=u,v=v):
  """ p = 2D surface in 3D space
  u,v = parameters on the surface
  """
  pu = diff(p,u)
  pv = diff(p,v)
  puu = diff(pu,u)
  puv = diff(pu,v)
  pvv = diff(pv,v)
  n = pu.cross(pv).normalize()
  E = pu.dot(pu) # 1st fundamental form
  F = pu.dot(pv)
  G = pv.dot(pv)
  L = puu.dot(n) # 2nd fundamental form
  M = puv.dot(n)
  N = pvv.dot(n)
  D = E*G - F**2
  # Weingarten's formula for
  K = (L*N - M**2)/D # Gauss curvature and
  H = (G*L - 2*F*M + E*N)/(2*D) # mean curvature
  return simplify(K), simplify(H)
