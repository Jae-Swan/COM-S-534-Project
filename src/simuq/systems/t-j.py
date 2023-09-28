__author__ = "Jae Swanepoel"

from simuq.environment import Fermion
from simuq.qsystem import QSystem
import numpy as np

# Pauli operators
# Not sure if these operators are defined elsewhere
pX = [[0, 1], [1, 0]]
pY = [[0, -1j], [1j, 0]]
pZ = [[1, 0], [0, -1]]
pops = [pX, pY, pZ]

L = 5  # number of lattice sites
N = 2 * L  # number of fermionic systems

qs = QSystem()
fermions = [Fermion(qs) for _ in range(N)]

sops = []  # spin operators
for k in range(L):
    u, d = 2 * k, 2 * k + 1
    sops[k] = [0.5 * (fermions[u].c * fermions[d].a + fermions[d].c * fermions[u].a),
               0.5j * (fermions[d].c * fermions[u].a - fermions[u].c * fermions[d].a),
               0.5 * (fermions[u].c * fermions[u].a - fermions[d].c * fermions[d].a)]

t = 0.1  # hopping amplitude
U = 0.3  # strength of on-site interaction
J = (t * t) / U  # exchange interaction
ht = 0  # t term
hJ = 0  # J term

# TODO I think the t term calculation may be wrong
for i in range(L - 1):
    for s in range(2):
        f1, f2 = 2 * i + s, 2 * (i + 1) + s
        ht += -t * (fermions[f1].c * fermions[f2].a + fermions[f2].c * fermions[f1].a)

for i in range (L):
    hJ += J * (sops[i] * sops[i + 1])

h = ht + hJ
qs.add_evolution(h, t)
