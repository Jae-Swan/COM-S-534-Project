from simuq import QSystem, Qubit
from simuq.systems import tJ
from simuq.systems import ising
from simuq.qutip import QuTiPProvider
from simuq.transformation import jw_transform

qtpp = QuTiPProvider()


def test_tJ():
    N = 3
    time = 0.1
    t = 0.1
    U = 0.3
    bin_len = N * 2

    # altering time
    print("testing changes over time")
    for i in range(1, 11):
        tJ_system = tJ.GentJ(N, time * i, t * i, U)
        print("\niteration ", str(i))
        qtpp.compile(tJ_system)
        qtpp.run()
        results = qtpp.results()
        for j in range(2 ** N):
            print(format(j, f'0{bin_len}b'), " : ", results[format(j, f'0{bin_len}b')])

    print("\ntesting changes over t")
    # testing altering t
    for i in range(1, 11):
        tJ_system = tJ.GentJ(N, time, t * i, U)
        print("\niteration ", str(i))
        qtpp.compile(tJ_system)
        qtpp.run()
        results = qtpp.results()
        for j in range(2 ** N):
            print(format(j, f'0{bin_len}b'), " : ", results[format(j, f'0{bin_len}b')])

    print("\ntesting changes over U")
    # altering U
    for i in range(1, 11):
        tJ_system = tJ.GentJ(N, time, t, U * i)
        print("\niteration ", str(i))
        qtpp.compile(tJ_system)
        qtpp.run()
        results = qtpp.results()
        for j in range(2 ** N):
            print(format(j, f'0{bin_len}b'), " : ", results[format(j, f'0{bin_len}b')])


def test_ising():
    N = 2
    time = 0.1
    J = 0.1
    h = 0.1
    bin_len = N

    print("testing changes over time:")
    for i in range(1, 11):
        ising_system = ising.GenQS(N, time * i, J, h)
        print("\niteration ", str(i))
        qtpp.compile(ising_system)
        qtpp.run()
        results = qtpp.results()
        for j in range(2 ** N):
            print(format(j, f'0{bin_len}b'), " : ", results[format(j, f'0{bin_len}b')])

    print("\ntesting changes over J:")
    for i in range(1, 11):
        ising_system = ising.GenQS(N, time, J * i, h)
        print("\niteration ", str(i))
        qtpp.compile(ising_system)
        qtpp.run()
        results = qtpp.results()
        for j in range(2 ** N):
            print(format(j, f'0{bin_len}b'), " : ", results[format(j, f'0{bin_len}b')])

    print("\ntesting changes over h:")
    for i in range(1, 11):
        ising_system = ising.GenQS(N, time, J, h * i)
        print("\niteration ", str(i))
        qtpp.compile(ising_system)
        qtpp.run()
        results = qtpp.results()
        for j in range(2 ** N):
            print(format(j, f'0{bin_len}b'), " : ", results[format(j, f'0{bin_len}b')])


def test_jw_transform():
    qs = tJ.GentJ(2, 0.1, 0.1, 0.3)
    print("Fermionic Hamiltonian:")
    for i in range(len(qs.evos[0][0].ham)):
        print(qs.evos[0][0].ham[i])
    new_qs, new_sites = jw_transform(qs)
    print("sites: ", new_qs.print_sites())
    print("\nQubit Hamiltonian:")
    for i in range(len(new_qs.evos[0][0].ham)):
        print(new_qs.evos[0][0].ham[i])
    qtpp.compile(new_qs)
    qtpp.run()
    result = qtpp.results()
    print("\nSite Frequencies: ")
    for j in range(2 ** 4):
        print(format(j, f'0{4}b'), " : ", result[format(j, f'0{4}b')])


def test_compare_tJ_ising():
    print("t-J System:")
    test_jw_transform()
    print("\nIsing System:")
    print("Qubit Hamiltonian")
    ising_system = ising.GenQS(4, 0.1, 0.1, 0.1)
    for i in range(len(ising_system.evos[0][0].ham)):
        print(ising_system.evos[0][0].ham[i])
    qtpp.compile(ising_system)
    qtpp.run()
    results = qtpp.results()
    print("\nSite Frequencies: ")
    for i in range(2 ** 4):
        print(format(i, f'0{4}b'), " : ", results[format(i, f'0{4}b')])
