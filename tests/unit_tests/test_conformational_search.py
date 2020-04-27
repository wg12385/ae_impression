
import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')

from autoenrich.conformational_search import conformational_search as csearch

import numpy as np

# Run torsional search for ethanol
def test_torsional_search():

    smiles = 'CCO'
    iterations = 2

    xyzs, energies = csearch.torsional_search(smiles, iterations)

    # xyz array contains [x, y, z] entries
    assert len(xyzs[0][0]) == 3
    # All xyz arrays are the same length
    for i in range(len(xyzs)):
        for j in range(i, len(xyzs)):
            assert len(xyzs[i]) == len(xyzs[j])
    # At least some energies are non-zero
    assert np.sum(energies) != 0


def test_select_conformers(size=10):

    xyzs = []
    energies = []
    for i in range(size):
        xyz = []
        for t in range(5):
            xyz.append([np.random.rand(), np.random.rand(), np.random.rand()])
        xyzs.append(xyz)
        energies.append(np.random.rand()*999)

    xyzs, energies = csearch.select_conformers(xyzs, energies, maxconfs=3)
    # xyz array contains [x, y, z] entries
    assert len(xyzs[0][0]) == 3
    # only 3 conformers remain, specified by maxconfs
    assert len(xyzs) == 3
    # All xyz arrays are the same length
    for i in range(len(xyzs)):
        for j in range(i, len(xyzs)):
            assert len(xyzs[i]) == len(xyzs[j])
    # At least some energies are non-zero
    assert np.sum(energies) != 0

    xyzs = []
    energies = []
    for i in range(size):
        xyz = []
        for t in range(5):
            xyz.append([np.random.rand(), np.random.rand(), np.random.rand()])
        xyzs.append(xyz)
        energies.append(np.random.rand()*999)

    xyzs, energies = csearch.select_conformers(xyzs, energies, Ethresh=sorted(energies)[2])
    # xyz array contains [x, y, z] entries
    assert len(xyzs[0][0]) == 3
    # only 3 conformers remain, specified by selecting 3rd lowest energy as Ethresh
    assert len(xyzs) == 3
    # All xyz arrays are the same length
    for i in range(len(xyzs)):
        for j in range(i, len(xyzs)):
            assert len(xyzs[i]) == len(xyzs[j])
    # At least some energies are non-zero
    assert np.sum(energies) != 0














###
