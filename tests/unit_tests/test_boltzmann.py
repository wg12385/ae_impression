import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')

from autoenrich.boltzmann.averaging import boltzmann_shift, boltzmann_coupling
from autoenrich.boltzmann.population import get_pop_array
from test_generators.dummy_mol import get_random_mol_with_confs

import numpy as np

def test_boltzmann_shift():

    mol = get_random_mol_with_confs()

    shift_array, shift_var = boltzmann_shift(mol.conformers)
    # averaged array is the same shape as molecule array
    assert shift_array.shape == mol.shift.shape
    # Shift array is correct shape
    assert shift_array.shape[0] == len(mol.types)
    # Shift variance array is correct shape
    assert shift_var.shape == shift_array.shape
    # At least some shifts are non-zero
    assert np.sum(shift_array) != 0


def test_boltzmann_coupling():

    mol = get_random_mol_with_confs()

    coupling_array, coupling_var = boltzmann_coupling(mol.conformers)
    # averaged array is same shape as molecule array
    assert coupling_array.shape == mol.coupling.shape
    # Coupling array is correct shape
    assert coupling_array.shape[0] == len(mol.types)
    assert coupling_array.shape[1] == len(mol.types)
    # Coupling variance array is same shape as coupling array
    assert coupling_var.shape == coupling_array.shape
    # At least some couplings are non-zero
    assert np.sum(coupling_array) != 0


def test_get_pop_array():

    mol = get_random_mol_with_confs()

    pop_array = get_pop_array(mol.conformers)
    # Population roughly sums to 1
    assert np.absolute(1 - np.sum(pop_array)) < 0.001
    # Population array contains an entry for each conformer
    assert pop_array.shape[0] == len(mol.conformers)
