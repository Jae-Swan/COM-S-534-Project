__author__ = "Jae Swanepoel"

from simuq.environment import Fermion
from simuq.qsystem import QSystem


def GentJ(num_sites):
    L = num_sites  # number of lattice sites
    N = 2 * L  # number of fermionic systems

    qs = QSystem()
    fermions = [Fermion(qs) for _ in range(N)]

    sops = [[0, 0, 0] for _ in range(L)]  # spin operators
    for i in range(L):
        u, d = fermions[i], fermions[i + 1]
        sops[i] = [0.5 * (u.c * d.a + d.c * u.a),
                   0.5j * (d.c * u.a - u.c * d.a),
                   0.5 * (u.c * u.a - d.c * d.a)]

    nops = [0 for _ in range(L)]
    for i in range(L):
        nops[i] += fermions[i * 2].c * fermions[i * 2].a
        nops[i] += fermions[i * 2 + 1].c * fermions[i * 2 + 1].a

    t = 0.1  # hopping amplitude
    U = 0.3  # strength of on-site interaction
    J = 2 * (t * t) / U  # exchange interaction
    h = 0  # Hamiltonian

    for i in range(L - 1):
        iu, id, ju, jd = fermions[2 * i], fermions[2 * i + 1], fermions[2 * (i + 1)], fermions[2 * (i + 1) + 1]
        h += (-1 * t) * ((iu.c * ju.a) + (id.c * jd.a) + (ju.c * iu.a) + (jd.c * id.a))  # t term

        Si, Sj = sops[i], sops[i + 1]
        ni, nj = nops[i], nops[i + 1]
        h += J * ((Si[0] * Sj[0]) + (Si[1] * Sj[1]) + (Si[2] * Sj[2]) - (0.25 * ni * nj))  # J term

    qs.add_evolution(h, t)
    return qs
