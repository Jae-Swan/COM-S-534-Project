__author__ = "Jae Swanepoel"

from simuq.environment import Fermion
from simuq.qsystem import QSystem


def GentJ(num_sites, time, _t, _U):
    L = num_sites  # number of lattice sites
    N = 2 * L  # number of fermionic systems

    qs = QSystem()
    fermions = [Fermion(qs) for _ in range(N)]

    # spin operators
    sops = [[0, 0, 0] for _ in range(L)]
    for i in range(L):
        u, d = fermions[i], fermions[i + 1]
        sops[i] = [((u.c * d.a) + (d.c * u.a)) / 2,
                   1j * ((d.c * u.a) - (u.c * d.a)) / 2,
                   ((u.c * u.a) - (d.c * d.a)) / 2]

    # number operator
    nops = [0 for _ in range(L)]
    for i in range(L):
        nops[i] = fermions[i * 2].c * fermions[i * 2].a
        nops[i] += fermions[i * 2 + 1].c * fermions[i * 2 + 1].a

    t = _t  # hopping amplitude
    U = _U  # strength of on-site interaction
    J = 2 * (t**2) / U  # exchange interaction
    ht = 0  # t-term
    hJ = 0  # J-term

    for i in range(L - 1):
        iu, _id, ju, jd = fermions[2 * i], fermions[(2 * i) + 1], fermions[2 * (i + 1)], fermions[(2 * (i + 1)) + 1]
        ht += (iu.c * ju.a) + (_id.c * jd.a) + (ju.c * iu.a) + (jd.c * _id.a)  # t term

        Si, Sj = sops[i], sops[i + 1]
        ni, nj = nops[i], nops[i + 1]
        hJ += ((Si[0] * Sj[0]) + (Si[1] * Sj[1]) + (Si[2] * Sj[2])) - (ni * nj / 4)  # J term

    h = (-t * ht) + (J * hJ)
    qs.add_evolution(h, time)
    return qs
