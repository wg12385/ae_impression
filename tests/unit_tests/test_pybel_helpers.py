


import pybel as pyb
import numpy as np

from test_generators.dummy_mol import get_ethane_mol

from autoenrich.pybel_helpers import pybel_bonds as pbonds
from autoenrich.pybel_helpers import pybel_read as pread
from autoenrich.pybel_helpers import pybel_split as psplit


def test_mol_find_all_paths():

    file = 'tests/test_store/ethane.xyz'
    pmol = next(pyb.readfile('xyz', file))
    mol = get_ethane_mol()
    size = len(mol.types)

    for i in range(size):
        for j in range(size):
            for k in range(size):
                paths = pbonds.mol_find_all_paths(pmol, i, j, k)

                if i == j:
                    assert len(paths) == 0
                else:
                    for path in paths:
                        if len(path) == 0:
                            continue
                        assert path[0] == i
                        assert path[-1] == j
                        assert len(set(path)) == len(path)


                    if mol.coupling_len[i][j] == k:
                        assert len(paths[0]) == k + 1

def test_mol_get_bond_table():

    file = 'tests/test_store/ethane.xyz'
    pmol = next(pyb.readfile('xyz', file))
    mol = get_ethane_mol()
    size = len(mol.types)

    bond_table = pbonds.mol_get_bond_table(pmol)
    assert np.array_equal(bond_table, mol.conn)


def test_get_coupling_lengths():

    file = 'tests/test_store/ethane.xyz'
    pmol = next(pyb.readfile('xyz', file))
    mol = get_ethane_mol()
    size = len(mol.types)

    coupling_len = pbonds.get_coupling_lengths(pmol, mol.types, maxlen=4)
    assert np.array_equal(coupling_len, mol.coupling_len)


def test_mol_read_type():

    file = 'tests/test_store/ethane.xyz'
    pmol = next(pyb.readfile('xyz', file))
    mol = get_ethane_mol()
    size = len(mol.types)

    type_list, type_array = pread.mol_read_type(pmol)
    assert np.array_equal(mol.types, type_array)
    assert np.array_equal(np.asarray([x for x, type in enumerate(type_list) if "H" == type]),
                            np.where(type_array==1)[0])
    assert np.array_equal(np.asarray([x for x, type in enumerate(type_list) if "C" == type]),
                            np.where(type_array==6)[0])

def test_mol_getatoms():

    file = 'tests/test_store/ethane.xyz'
    pmol = next(pyb.readfile('xyz', file))
    mol = get_ethane_mol()
    size = len(mol.types)
    assert size == pread.mol_getatoms(pmol)

def test_mol_read_xyz():

    file = 'tests/test_store/ethane.xyz'
    pmol = next(pyb.readfile('xyz', file))
    mol = get_ethane_mol()

    assert np.array_equal(mol.xyz, pread.mol_read_xyz(pmol))

def test_mol_read_dist():

    file = 'tests/test_store/ethane.xyz'
    pmol = next(pyb.readfile('xyz', file))
    mol = get_ethane_mol()
    size = len(mol.types)
    d_array = pread.mol_read_dist(pmol)
    for i in range(size):
        for j in range(size):
            assert d_array[i][j] == np.absolute(np.linalg.norm(mol.xyz[i] - mol.xyz[j]))


def test_mol_iswhole():

    file = 'tests/test_store/ethane.xyz'
    pmol = next(pyb.readfile('xyz', file))
    mol = get_ethane_mol()
    size = len(mol.types)

    assert psplit.mol_iswhole(pmol)


def test_mol_splitmol():

    file = 'tests/test_store/ethane.xyz'
    pmol = next(pyb.readfile('xyz', file))
    mol = get_ethane_mol()
    size = len(mol.types)

    print(psplit.mol_splitmol(pmol))

    assert len(psplit.mol_splitmol(pmol)[0]) == size
    assert len(psplit.mol_splitmol(pmol)[1]) == 0
    assert len(psplit.mol_splitmol(pmol)) == 2

    # doesnt test splitting a mol,
    # just tests that trying to split a whole molecule does nothing




#
