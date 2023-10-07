__author__ = "Jae Swanepoel"

from simuq.environment import Fermion
from simuq.qsystem import QSystem
from simuq.qutip import QuTiPProvider

qtpp = QuTiPProvider()

L = 5  # number of lattice sites
N = 2 * L  # number of fermionic systems

qs = QSystem()
fermions = [Fermion(qs) for _ in range(N)]

sops = [[0, 0, 0] for _ in range(L)]  # spin operators
for i in range(L):
    u, d = fermions[i], fermions[i + 1]
    sops[i] = [0.5 * (u.c * d.a + d.c * u.a),
               0.5j * (d.c * u.a - u.c * d.a),
               0.5 * (u.c * u.a - d.c * d.a)]

particle_numbers = [0 for _ in range(L)]
for i in range(L):
    particle_numbers[i] = fermions[i * 2].c * fermions[i * 2].a + fermions[i * 2 + 1].c * fermions[i * 2 + 1].a

t = 0.1  # hopping amplitude
U = 0.3  # strength of on-site interaction
J = (t * t) / U  # exchange interaction
h = 0  # Hamiltonian

for i in range(L - 1):
    iu, id, ju, jd = fermions[2 * i], fermions[2 * i + 1], fermions[2 * (i + 1)], fermions[2 * (i + 1) + 1]
    h += (-1 * t) * (iu.c * ju.a + id.c * jd.a + ju.c * iu.a + jd.c * id.a)  # t term
    Si, Sj = sops[i], sops[i + 1]
    h += J * ((Si[0] * Sj[0]) + (Si[1] * Sj[1]) + (Si[2] * Sj[2]) - (0.25 * particle_numbers[i] * particle_numbers[i + 1]))  # J term = Si.Sj

qs.add_evolution(h, t)

qtpp.compile(qs)
qtpp.run()
res_cycle_gt = qtpp.results()
print(res_cycle_gt["0000000000"])
