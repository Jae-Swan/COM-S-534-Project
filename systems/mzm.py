import numpy as np
from simuq.qsystem import QSystem
from simuq.environment import qubit
from simuq.hamiltonian import Empty
# Minimal example of the braiding of Majorana zero mode
qs = QSystem()
ql = [qubit(qs) for i in range(3)]

Jmax = 1.0
alpha = Jmax
tau = 3.3 / Jmax

hbase = Empty
for q_id in range(3):
    hbase += alpha * ql[q_id].Z

h01 = hbase + Jmax * ql[0].Y * ql[1].X
h12 = hbase + Jmax * ql[1].Y * ql[2].X
h20 = hbase + Jmax * ql[0].Y * ql[1].Z * ql[2].X


## qs with 3-qubit pauli gates
#for step_id in range(2): # repeat twice
#    qs.add_evolution(h12, tau)
#    qs.add_evolution(h20, tau)
#    qs.add_evolution(h01, tau)

happrox0 = ql[0].Y * ql[1].X
happrox1 = ql[1].Y * ql[2].X



# Using 1st order trotter formula and approximation in Eq (11), arxiv:2003.06886
for step_id in range(2): # repeat twice
    qs.add_evolution(h12, tau)
    qs.add_evolution(hbase, tau)
    qs.add_evolution(-happrox0, np.sqrt(tau*Jmax/2))
    qs.add_evolution(happrox1, np.sqrt(tau*Jmax/2))
    qs.add_evolution(happrox0, np.sqrt(tau*Jmax/2))
    qs.add_evolution(-happrox1, np.sqrt(tau*Jmax/2))
    qs.add_evolution(h01, tau)



