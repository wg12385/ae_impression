import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../../')

from autoenrich.molecule.nmrmol import nmrmol
from autoenrich.molecule.molecule import molecule
from autoenrich.molecule.conformer import conformer
import numpy as np


def get_random_mol(size=10):

    mol = nmrmol(0)

    mol.types = np.random.choice([1, 6, 7, 8], size=size)
    mol.xyz = np.random.rand(size, 3)
    conn = np.zeros((size, size), dtype=np.int32)
    for x in range(size):
        for y in range(x, size):
            conn[x][y] = np.random.choice([0, 1, 2])
            conn[y][x] = conn[x][y]
    mol.conn = conn

    mol.shift = np.random.rand(size)
    mol.shift_var = np.random.rand(size)
    mol.coupling = np.random.rand(size, size)
    mol.coupling_var = np.random.rand(size, size)

    mol.coupling_len = mol.conn


    return mol


def get_random_mol_with_confs(size=5):

    mol = molecule(0)

    mol.types = np.random.choice([1, 6, 7, 8], size=size)
    mol.xyz = np.random.rand(size, 3)
    conn = np.zeros((size, size), dtype=np.int32)
    for x in range(size):
        for y in range(x, size):
            conn[x][y] = np.random.choice([0, 1, 2])
            conn[y][x] = conn[x][y]
    mol.conn = conn

    mol.shift = np.random.rand(size)
    mol.shift_var = np.random.rand(size)
    mol.coupling = np.random.rand(size, size)
    mol.coupling_var = np.random.rand(size, size)

    mol.coupling_len = []

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
