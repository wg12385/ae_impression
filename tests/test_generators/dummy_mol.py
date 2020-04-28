import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../../')

from autoenrich.molecule.nmrmol import nmrmol
from autoenrich.molecule.molecule import molecule
from autoenrich.molecule.conformer import conformer
import numpy as np

def mol_is_ethane(mol):
    '''
    assert np.array_equal(mol.xyz, [[-0.4125,  0.0000,  0.0000],
                     [ 0.4125,  0.0000,  0.0000],
                     [-0.9475,  0.9266,  0.0000],
                     [-0.9475, -0.9266,  0.0000],
                     [-1.4825, -0.0000,  0.0000],
                     [ 0.9475, -0.9266,  0.0000],
                     [ 0.9475,  0.9266,  0.0000],
                     [ 1.4825,  0.0000,  0.0000]])
    '''
    assert np.array_equal(mol.types, [6, 6, 1, 1, 1, 1, 1, 1])

    assert np.array_equal(mol.conn, [[0, 1, 1, 1, 1, 0, 0, 0],
                                         [1, 0, 0, 0, 0, 1, 1, 1],
                                         [1, 0, 0, 0, 0, 0, 0, 0],
                                         [1, 0, 0, 0, 0, 0, 0, 0],
                                         [1, 0, 0, 0, 0, 0, 0, 0],
                                         [0, 1, 0, 0, 0, 0, 0, 0],
                                         [0, 1, 0, 0, 0, 0, 0, 0],
                                         [0, 1, 0, 0, 0, 0, 0, 0]])

    assert np.array_equal(mol.coupling_len, [[0, 1, 1, 1, 1, 2, 2, 2],
                                         [1, 0, 2, 2, 2, 1, 1, 1],
                                         [1, 2, 0, 2, 2, 3, 3, 3],
                                         [1, 2, 2, 0, 2, 3, 3, 3],
                                         [1, 2, 2, 2, 0, 3, 3, 3],
                                         [2, 1, 3, 3, 3, 0, 2, 2],
                                         [2, 1, 3, 3, 3, 2, 0, 2],
                                         [2, 1, 3, 3, 3, 2, 2, 0]])

    return True

def get_ethane_mol():

    mol = nmrmol(molid='ethane')

    mol.types = np.asarray([6, 6, 1, 1, 1, 1, 1, 1])
    size = len(mol.types)
    mol.xyz = np.asarray([[-0.4125,  0.0000,  0.0000],
                     [ 0.4125,  0.0000,  0.0000],
                     [-0.9475,  0.9266,  0.0000],
                     [-0.9475, -0.9266,  0.0000],
                     [-1.4825, -0.0000,  0.0000],
                     [ 0.9475, -0.9266,  0.0000],
                     [ 0.9475,  0.9266,  0.0000],
                     [ 1.4825,  0.0000,  0.0000]])

    mol.conn = np.asarray([[0, 1, 1, 1, 1, 0, 0, 0],
                         [1, 0, 0, 0, 0, 1, 1, 1],
                         [1, 0, 0, 0, 0, 0, 0, 0],
                         [1, 0, 0, 0, 0, 0, 0, 0],
                         [1, 0, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 0, 0]])

    mol.coupling_len = np.asarray([[0, 1, 1, 1, 1, 2, 2, 2],
                                     [1, 0, 2, 2, 2, 1, 1, 1],
                                     [1, 2, 0, 2, 2, 3, 3, 3],
                                     [1, 2, 2, 0, 2, 3, 3, 3],
                                     [1, 2, 2, 2, 0, 3, 3, 3],
                                     [2, 1, 3, 3, 3, 0, 2, 2],
                                     [2, 1, 3, 3, 3, 2, 0, 2],
                                     [2, 1, 3, 3, 3, 2, 2, 0]])

    mol.shift = np.asarray([223.81,  223.81,   32.002,  32.002,  32.005,  32.002, 32.005,  32.002])
    mol.shift_var = np.zeros(size, dtype=np.float64)
    mol.coupling = np.asarray([[ 0.0000,  38.796, 114.462, 115.871, 113.801, -4.644 ,  -4.576,  -4.558],
                                         [ 38.796,   0.000,  -4.558,  -4.644,  -4.576, 115.87 , 113.801, 114.462],
                                         [114.462,  -4.558,   0.000, -13.961, -13.644,   3.174,   3.161,  12.579],
                                         [115.871,  -4.644, -13.961,   0.000, -13.906,  12.673,   3.181,   3.182],
                                         [113.801,  -4.576, -13.644, -13.906,   0.000,   3.178,  12.524,   3.167],
                                         [ -4.644, 115.87 ,   3.174,  12.673,   3.178,   0.000, -13.906, -13.971],
                                         [ -4.576, 113.801,   3.161,   3.181,  12.524, -13.906,   0.000, -13.658],
                                         [ -4.558, 114.462,  12.579,   3.182,   3.167, -13.971, -13.658,   0.000]])
    mol.coupling_var = np.zeros((size, size), dtype=np.float64)

    return mol



def get_random_ethane():
    # Gets an ethane molecule with slightly varied xyz and NMR parameters
    mol = get_ethane_mol()

    size = len(mol.types)
    mol.xyz = mol.xyz * (np.random.rand(size, 3)*0.02+0.99)
    mol.shift = mol.shift * (np.random.rand(size)*0.02+0.99)
    mol.shift_var = np.random.rand(size)
    mol.coupling = mol.coupling * (np.random.rand(size, size)*0.02+0.99)
    mol.coupling_var = np.random.rand(size, size)


    return mol

def get_random_mol(size=10):

    mol = molecule(0)
    mol.types = np.random.choice([1, 6, 7, 8], size=size)
    mol.xyz = np.random.rand(size, 3)
    conn = np.zeros((size, size), dtype=np.int32)
    for x in range(size):
        for y in range(x, size):
            conn[x][y] = np.random.choice([0, 1, 2, 3, 4])
            conn[y][x] = conn[x][y]
    mol.conn = conn

    mol.shift = np.random.rand(size)
    mol.shift_var = np.random.rand(size)
    mol.coupling = np.random.rand(size, size)
    mol.coupling_var = np.random.rand(size, size)

    mol.coupling_len = conn

    return mol

def get_random_mol_with_confs(size=5):

    mol = get_random_mol(size=size)

    for i in range(size):
        conf = conformer(i)

        conf.types = mol.types
        conf.conn = mol.conn

        conf.opt_status = 'successful'
        conf.nmr_status = 'successful'
        conf.energy = np.random.rand() * 999
        conf.pop = np.random.rand()

        conf.xyz = np.random.rand(size, 3)
        conf.shift = np.random.rand(size)
        conf.shift_var = np.random.rand(size)
        conf.coupling = np.random.rand(size, size)
        conf.coupling_var = np.random.rand(size, size)

        mol.conformers.append(conf)


    return mol










#
