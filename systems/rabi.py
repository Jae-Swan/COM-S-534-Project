from simuq.qsystem import QSystem
from simuq.environment import qubit, fock
from simuq.qmachine import *

qs = QSystem()
q = qubit(qs)
n_step = 5
T = 1.
omega0 = 0.1
omega1 = 0.3
omega = 0.2
for step in range(n_step) :
    t = step * T / n_step
    h = -omega0 / 2 * q.Z() - omega1 / 2 * (math.cos(omega * t) * q.X() - math.sin(omega * t) * q.Y())
    qs.add_evolution(h, T / n_step)

mach = QMachine()
q0 = qubit(mach)
L0 = SignalLine(mach)
ins1 = Instruction(L0, 'native', 'L0_XY')
a = LocalVar(ins1)
b = LocalVar(ins1)
ins1.set_ham(b * (Expression.cos(a) * q0.X() + Expression.sin(a) * q0.Y()))
ins2 = Instruction(L0, 'native', 'L0_Z')
b = LocalVar(ins2)
ins2.set_ham(b * q0.Z())