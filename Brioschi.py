from sympy import diff,simplify
from sympy.abc import u,v
from sympy.matrices import Matrix,det

def metric(p,u=u,v=v):
    """ p = 2D surface in 3D space
    u,v = parameters on the surface
    """
    pu = diff(p,u)
    pv = diff(p,v)
    E = pu.dot(pu)
    F = pu.dot(pv)
    G = pv.dot(pv)
    return E,F,G # 1st fundamental form

def Brioschi(metric,u=u,v=v):
    E,F,G = metric # 1st fundamental form 
    Eu = diff(E,u)
    Ev = diff(E,v)
    Fu = diff(F,u)
    Fv = diff(F,v)
    Gu = diff(G,u)
    Gv = diff(G,v)
    Evv = diff(Ev,v)
    Fuv = diff(Fu,v)
    Guu = diff(Gu,u)
    A = Matrix([[0, Eu/2, Fu - Ev/2],
                [Fv - Gu/2, E, F],
                [Gv/2, F, G]])
    B = Matrix([[Evv/2 - Fuv + Guu/2, Ev/2, Gu/2],
                [Ev/2, E, F],
                [Gu/2, F, G]])
    K = (det(A) - det(B))/(E*G - F**2)**2
    return simplify(K) # Gauss curvature
