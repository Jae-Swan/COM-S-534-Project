from simuq.systems import tJ
from simuq.systems import ising
from simuq.qutip import QuTiPProvider

qtpp = QuTiPProvider()

tJ_system = tJ.GentJ(5)
qtpp.compile(tJ_system)
qtpp.run()
res_cycle_gt = qtpp.results()
for i in range(len(res_cycle_gt)):
    print(format(i, '010b'), " : ", res_cycle_gt[format(i, '010b')])
