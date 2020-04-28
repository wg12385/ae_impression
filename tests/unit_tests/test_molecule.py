

import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')

import numpy as np

from autoenrich.molecule.nmrmol import nmrmol
from autoenrich.molecule.molecule import molecule
from autoenrich.molecule.conformer import conformer
from autoenrich.molecule.dataset import dataset

from test_generators.dummy_mol import mol_is_ethane

def test_nmrmol():

    mol = nmrmol(molid='ethane')

    for attr in [mol.types, mol.xyz, mol.conn, mol.shift,
                    mol.shift_var, mol.coupling, mol.coupling_var,
                        mol.coupling_len]:
        assert attr == []

    assert mol.path == ''
    assert mol.molid == 'ethane'
    assert mol.energy == -404.404

def test_molecule():

    mol = molecule(molid='ethane')

    for attr in [mol.types, mol.xyz, mol.conn, mol.shift,
                    mol.shift_var, mol.coupling, mol.coupling_var,
                        mol.coupling_len]:
        assert attr == []

    assert mol.path == ''
    assert mol.molid == 'ethane'
    assert mol.energy == -404.404

    assert mol.conformers == []
    assert mol.stage == 'init'

def test_conformer():

    mol = conformer(molid='ethane')

    for attr in [mol.types, mol.xyz, mol.conn, mol.shift,
                    mol.shift_var, mol.coupling, mol.coupling_var,
                        mol.coupling_len]:
        assert attr == []

    assert mol.path == ''
    assert mol.molid == 'ethane'
    assert mol.energy == 404.404

    assert mol.xyz_file == 'None'
    assert mol.opt_in == 'None'
    assert mol.opt_log == 'None'
    assert mol.opt_status == 'None'
    assert mol.nmr_in == 'None'
    assert mol.nmr_log == 'None'
    assert mol.nmr_status == 'None'
    assert mol.pop == 404.404
    assert mol.redundant == False

def test_nmrmol_read_structure():

    mol = nmrmol(molid='ethane')

    mol.read_structure('tests/test_store/ethane.xyz', 'xyz')

    assert mol_is_ethane(mol)

def test_nmrmol_read_opt():

    mol = nmrmol(molid='ethane')

    mol.read_opt('tests/test_store/TST_pass_orca_opt.log', 'orca')
    assert mol.energy == -1072.556525704736

    mol.read_opt('tests/test_store/TST_pass_g09_opt.log', 'g09')
    assert mol.energy == -1419.091348

    #mol.read_opt('tests/test_store/TST_pass_g16_opt.log', 'orca')
    #assert mol.energy == 2342.2

def test_nmrmol_read_nmr():

    mol = nmrmol(molid='ethane')

    mol.read_nmr('tests/test_store/ethane_orca_nmr.log', 'orca')

    assert np.array_equal(mol.shift, [223.81,  223.81,   32.002,  32.002,  32.005,  32.002, 32.005,  32.002])
    assert np.array_equal(mol.coupling, [[ 0.0000,  38.796, 114.462, 115.871, 113.801, -4.644 ,  -4.576,  -4.558],
                                         [ 38.796,   0.000,  -4.558,  -4.644,  -4.576, 115.87 , 113.801, 114.462],
                                         [114.462,  -4.558,   0.000, -13.961, -13.644,   3.174,   3.161,  12.579],
                                         [115.871,  -4.644, -13.961,   0.000, -13.906,  12.673,   3.181,   3.182],
                                         [113.801,  -4.576, -13.644, -13.906,   0.000,   3.178,  12.524,   3.167],
                                         [ -4.644, 115.87 ,   3.174,  12.673,   3.178,   0.000, -13.906, -13.971],
                                         [ -4.576, 113.801,   3.161,   3.181,  12.524, -13.906,   0.000, -13.658],
                                         [ -4.558, 114.462,  12.579,   3.182,   3.167, -13.971, -13.658,   0.000]])

    '''
    mol.read_nmr('tests/test_store/ethane_g09_nmr.log', 'g09')

    assert mol_is_ethane(mol)
    assert np.array_equal(mol.shift, [223.81,  223.81,   32.002,  32.002,  32.005,  32.002, 32.005,  32.002])
    assert np.array_equal(mol.coupling, [[ 0.0000,  38.796, 114.462, 115.871, 113.801, -4.644 ,  -4.576,  -4.558],
                                         [ 38.796,   0.000,  -4.558,  -4.644,  -4.576, 115.87 , 113.801, 114.462],
                                         [114.462,  -4.558,   0.000, -13.961, -13.644,   3.174,   3.161,  12.579],
                                         [115.871,  -4.644, -13.961,   0.000, -13.906,  12.673,   3.181,   3.182],
                                         [113.801,  -4.576, -13.644, -13.906,   0.000,   3.178,  12.524,   3.167],
                                         [ -4.644, 115.87 ,   3.174,  12.673,   3.178,   0.000, -13.906, -13.971],
                                         [ -4.576, 113.801,   3.161,   3.181,  12.524, -13.906,   0.000, -13.658],
                                         [ -4.558, 114.462,  12.579,   3.182,   3.167, -13.971, -13.658,   0.000]])

    mol.read_nmr('tests/test_store/ethane_g16_nmr.log', 'g09')

    assert mol_is_ethane(mol)
    assert np.array_equal(mol.shift, [223.81,  223.81,   32.002,  32.002,  32.005,  32.002, 32.005,  32.002])
    assert np.array_equal(mol.coupling, [[ 0.0000,  38.796, 114.462, 115.871, 113.801, -4.644 ,  -4.576,  -4.558],
                                         [ 38.796,   0.000,  -4.558,  -4.644,  -4.576, 115.87 , 113.801, 114.462],
                                         [114.462,  -4.558,   0.000, -13.961, -13.644,   3.174,   3.161,  12.579],
                                         [115.871,  -4.644, -13.961,   0.000, -13.906,  12.673,   3.181,   3.182],
                                         [113.801,  -4.576, -13.644, -13.906,   0.000,   3.178,  12.524,   3.167],
                                         [ -4.644, 115.87 ,   3.174,  12.673,   3.178,   0.000, -13.906, -13.971],
                                         [ -4.576, 113.801,   3.161,   3.181,  12.524, -13.906,   0.000, -13.658],
                                         [ -4.558, 114.462,  12.579,   3.182,   3.167, -13.971, -13.658,   0.000]])
    '''

    mol.read_nmr('tests/test_store/ethane.nmredata.sdf', 'nmredata')

    assert np.allclose(mol.shift, [223.849, 223.849, 32.006, 32.003,  32.012, 32.006, 32.004,  32.011], 0.001)
    assert np.allclose(mol.coupling, [[  0.000,  38.849, 114.872, 115.641, 114.354,  -4.636,  -4.627,  -4.609],
                                         [ 38.849,   0.000,  -4.607,  -4.65 ,  -4.615, 115.341, 114.951, 114.575],
                                         [114.872,  -4.607,   0.000, -14.014, -13.691,   9.479,   3.178,   6.319],
                                         [115.641,  -4.650, -14.014,   0.000, -13.953,   6.348,   9.489,   3.185],
                                         [114.354,  -4.615, -13.691, -13.953,   0.000,   3.187,   6.302,   9.464],
                                         [ -4.636, 115.341,   9.479,   6.348,   3.187,   0.000, -13.996, -13.800],
                                         [ -4.627, 114.951,   3.178,   9.489,   6.302, -13.996,   0.000, -13.870],
                                         [ -4.609, 114.575,   6.319,   3.185,   9.464, -13.800, -13.870,   0.000]], 0.001)

def test_scale_shifts():

    mol = nmrmol(molid='ethane')

    mol.read_structure('tests/test_store/ethane.xyz', 'xyz')
    mol.read_nmr('tests/test_store/ethane_orca_nmr.log', 'orca')

    mol.scale_shifts(scaling_factors={'H': [-1.1, 30.0], 'C': [-1.2, 30.2]})

    print(mol.shift)

    assert np.allclose(mol.shift, [-161.341, -161.341, -1.820, -1.820, -1.822, -1.820, -1.822, -1.820], 0.001)


def test_save_pickle():

    mol = nmrmol(molid='ethane')
    file = 'TMP_mol.pkl'
    mol.save_pickle(file)
    assert os.path.isfile(file)
    os.remove(file)













##
