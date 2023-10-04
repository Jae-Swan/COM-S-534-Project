__author__ = "Jae Swanepoel"

from simuq.environment import Fermion
from simuq.qsystem import QSystem
import numpy as np
from simuq.qutip import QuTiPProvider
from simuq.ionq import IonQProvider
iqp = IonQProvider("VEgePYu3LJtjGQBI4MFQrxmyfeGLRVBZ")

qpp = QuTiPProvider()


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

sops = [0 for _ in range(L)]  # spin operators
#TODO define sops

t = 0.1  # hopping amplitude
U = 0.3  # strength of on-site interaction
J = (t * t) / U  # exchange interaction
ht = 0  # t term
hJ = 0  # J term

for i in range(L - 1):
    iu, id, ju, jd = fermions[2 * i], fermions[2 * i + 1], fermions[2 * (i + 1)], fermions[2 * (i + 1) + 1]
    ht += -t * (iu.c * ju.a + id.c * jd.a + ju.c * iu.a + jd.c * id.a)

for i in range(L - 1):
    hJ += J * (sops[i] * sops[i + 1])

h = ht + hJ
qs.add_evolution(h, t)
print (qs)
qpp.compile(qs)
qpp.run()
qpp.results()
qpp.print_sites()

iqp.compile(qs)
iqp.run(on_simulator = True)
iqp.results()
iqp.print_sites()