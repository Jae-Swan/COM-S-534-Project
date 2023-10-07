from simuq.systems import tJ
from simuq.qutip import QuTiPProvider

qtpp = QuTiPProvider()

tJ_system = tJ.GentJ(5)
qtpp.compile(tJ_system)
qtpp.run()
res_cycle_gt = qtpp.results()
print(res_cycle_gt["0000000000"])