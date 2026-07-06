import numpy as np
import scipy.special as sp
from scipy.optimize import brentq

v = 2 * np.pi * 6e-6 * 0.0779 / 1065e-9
n2 = 1.4515
n1 = np.sqrt(0.0779**2-n2**2)

def neff(w, v, n1, n2):
    return n2*(1+(w/v)**2*(n1-n2)/(n1))

def u(l, m, v):
    u_inf = sp.jn_zeros(l, m)[-1]
    return u_inf*(v/(v+1))*(1-(u_inf**2)/(6*(v+1)**3)-(u_inf**4)/(20*(v+1)**5))

def psy(l, w):
    return sp.kn(l, w)**2/(sp.kn(l+1,w)*sp.kn(l-1,w))

def tau(u, v, w, l):
    return 1-(u**2/v**2)*(1-psy(l, w))

def w(v, u):
    return np.sqrt(v**2 - u**2)

def mode_bessel1(l, u):
    return (sp.jn(l, u)) / (u * sp.jn(l-1, u))

def mode_bessel1_2(l, w):
    return (sp.kv(l, w)) / (w * sp.kv(l-1, w))

def eq87(u, l, v):
    return mode_bessel1(l, u) + mode_bessel1_2(l, w(v, u))

for i, l in enumerate([0, 1, 2, 3, 4, 5]):
    for j, m in enumerate([1, 2, 3, 4, 5]):
        try:
            if 0<u(l, m, v)<v:
                print(f"u({l}, {m}, {v}) = {u(l, m, v)} = {neff(w(v, u(l, m, v)), v, n1, n2)} = {tau(u(l, m, v), v, w(v, u(l, m, v)), l)}")
        except:
            continue

def find_roots(l, v, u_min, u_max, num_points=100000, tol=1e-12):
    u_grid = np.linspace(u_min, u_max, num_points)
    eq_values = eq87(u_grid, l, v)
    roots = []
    for i in range(len(u_grid) - 1):
        if eq_values[i] * eq_values[i + 1] < 0:
            root = brentq(eq87, u_grid[i], u_grid[i + 1], args=(l, v), xtol=tol)
            roots.append(root)
    return roots #retourne des duos de racine (début de l'intervalle, intersection avec l'autre Bessel)

for l in range(0, 9):
    print (l)
    u_min = 0.001
    u_max = v
    roots = find_roots(l, v, u_min, u_max)
    print(f"Found {len(roots)} roots:")
    for i, root in enumerate(roots):
        #on ne considère que les racines paires, sauf pour l'itération l=0 qui n'en voit qu'une sur la première 
        print(f"Root {i + 1}: u = {root}, neff = {neff(w(v, root), v, n1, n2)}, tau = {tau(root, v, w(v, root), l)}")