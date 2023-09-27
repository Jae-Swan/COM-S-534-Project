__author__ = "Jae Swanepoel"

from simuq.environment import Fermion
from simuq.qsystem import QSystem
from scipy.constants import hbar

# Pauli operators
# Not sure if these operators are defined elsewhere
pX = [[0, 1], [1, 0]]
pY = [[0, -1j], [1j, 0]]
pZ = [[1, 0], [0, -1]]
pops = [pX, pY, pZ]

L = 5      # number of lattice sites
N = 2 * L  # number of fermionic systems

qs = QSystem()
fermions = [Fermion(qs) for _ in range(N)]

sops = []  # spin operators
for k in range(L):
    sops[k] = fermions[k].c * pops * fermions[k].a

t = 0.1          # hopping amplitude
U = 0.3          # strength of on-site interaction
J = (t * t) / U  # exchange interaction
hh = 0           # hopping term
ho = 0           # on-site interaction term

for i in range(L - 1):
    for s in range(2):
        f1, f2 = 2 * i + s, 2 * (i + 1) + s
        hh += (-t * (fermions[f1].c * fermions[f2].a + fermions[f2].c * fermions[f1].a)) + (J * (sops[f1] * sops[f2]))
#                                ^ t term ^                                                     ^ J term ^

for i in range(L - 1):
    ho += J * (sops[i] * sops[i + 1])

h = hh + ho
qs.add_evolution(h, t)
