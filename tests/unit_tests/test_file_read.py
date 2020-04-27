import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')

from autoenrich.file_read import g09_read as g09read
from autoenrich.file_read import g16_read as g16read
from autoenrich.file_read import nmredata_read as nmreread
from autoenrich.file_read import orca_read as orcaread
from autoenrich.file_read import structure_read as strucread


import numpy as np

def test_g09get_opt_status():

    print(sys.path)

    file = 'tests/test_store/TST_pass_g09_opt.log'
    status = g09read.get_opt_status(file)
    assert status == 'successful'

    file = 'tests/test_store/TST_fail_g09_opt.log'
    status = g09read.get_opt_status(file)
    assert status == 'failed'

def test_g09get_nmr_status():

    file = 'tests/test_store/TST_pass_g09_nmr.log'
    status = g09read.get_nmr_status(file)
    assert status == 'successful'

    file = 'tests/test_store/TST_fail_g09_nmr.log'
    status = g09read.get_nmr_status(file)
    assert status == 'unknown'

def test_g09read_opt():

    file = 'tests/test_store/TST_pass_g09_opt.log'
    energy = g09read.read_opt(file)
    assert energy == -1419.091348

def test_g09read_functional():

    file = 'tests/test_store/TST_pass_g09_nmr.log'
    functional, basisset = g09read.read_functional(file)

    assert functional == '(giao,spinspin,mixed)wb97xd'
    assert basisset == '6-311g(d,p)'

def test_g09read_nmr():

    file = 'tests/test_store/TST_pass_g09_nmr.log'

    shift, coupling = g09read.read_nmr(file)

    print(shift, coupling)



    shiftVL = np.asarray([101.9613, 140.5238, 106.7185, 143.8829, 162.7263, 148.0998, 111.091 , 117.9161,
                     124.7025, 120.0261, 112.2024,  18.8315,  42.4334, 260.3007, 834.0869, 257.5997,
                      29.2217, 169.589 , 248.6725,  160.975, -69.2941,  64.3367,  28.1209,  28.7921,
                      30.4229,  30.2125,  30.2792,  30.1523,  28.1316,  28.784 ,  31.0804,  30.3319,
                      29.3494,  30.6823,  30.5061,  30.8292,  30.5077,  26.3801,  25.4293,  28.5733])
    assert np.array_equal(shift, shiftVL)

    # Too many values to construct array so check a few
    assert coupling[0][0] == 0.0
    assert coupling[0][1] == 34.6201
    assert coupling[-1][-1] == 0.0
    assert coupling[4][34] == -0.19978

'''
def test_g16get_opt_status():

    file = 'tests/test_store/TST_pass_g16_opt.log'
    status = g16read.get_opt_status(file)
    assert status == 'successful'

    file = 'tests/test_store/TST_fail_g16_opt.log'
    status = g16read.get_opt_status(file)
    assert status == 'failed'

def test_g16get_nmr_status():

    file = 'tests/test_store/TST_pass_g16_nmr.log'
    status = g16read.get_nmr_status(file)
    assert status == 'successful'

    file = 'tests/test_store/TST_fail_g16_nmr.log'
    status = g16read.get_nmr_status(file)
    assert status == 'unknown'

def test_g16read_opt():

    file = 'tests/test_store/TST_pass_g16_opt.log'
    energy = g16read.read_opt(file)
    assert energy == -1419.091348

def test_g16read_functional():

    file = 'tests/test_store/TST_pass_g16_nmr.log'
    functional, basisset = g16read.read_functional(file)

    assert functional == '(giao,spinspin,mixed)wb97xd'
    assert basisset == '6-311g(d,p)'

def test_g16read_nmr():

    file = 'tests/test_store/TST_pass_g16_nmr.log'

    shift, coupling = g16read.read_nmr(file)

    print(shift, coupling)



    shiftVL = np.asarray([101.9613, 140.5238, 106.7185, 143.8829, 162.7263, 148.0998, 111.091 , 117.9161,
                     124.7025, 120.0261, 112.2024,  18.8315,  42.4334, 260.3007, 834.0869, 257.5997,
                      29.2217, 169.589 , 248.6725,  160.975, -69.2941,  64.3367,  28.1209,  28.7921,
                      30.4229,  30.2125,  30.2792,  30.1523,  28.1316,  28.784 ,  31.0804,  30.3319,
                      29.3494,  30.6823,  30.5061,  30.8292,  30.5077,  26.3801,  25.4293,  28.5733])
    assert np.array_equal(shift, shiftVL)

    # Too many values to construct array so check a few
    assert coupling[0][0] == 0.0
    assert coupling[0][1] == 34.6201
    assert coupling[-1][-1] == 0.0
    assert coupling[4][34] == -0.19978
'''


def test_nmreread_nmr():

    file = 'tests/test_store/TST.nmredata.sdf'

    shift_array, shift_var, coupling_array, coupling_var, coupling_len = nmreread.read_nmr(file)

    assert coupling_var.all() == 0
    assert shift_var.all() == 0

    shift = [143.48170367, 126.49607271, 138.00411311, 155.96715011, 126.40840368,
             168.29285238, 186.87592914, 183.36560189, 100.12294058, 263.12286246,
              97.97775459, 178.29262857, 137.85949488,  42.09028552, 236.50090923,
              37.15290164, 367.96346644, 184.95539942, 132.08494036, 221.66173856,
             219.92242374, 206.27279126,  39.41897739,  34.44230541, 214.78224299,
             -51.24261098,  50.9902238 , 110.43377212, 148.8344559 ,  38.35349795,
              34.02211866,  26.72205086,  31.53141487,  33.43826347,  29.57039314,
              32.57140199,  31.85189784,  24.69238035,  30.90415101,  31.7175425,
              37.16425071,  45.26465065,  33.60985415,  35.49155789,  42.79791116,
              32.40941756,  34.30290677]
    assert np.array_equal(shift_array, shift)

    coupling0 =  np.asarray([0.00000000e+00,  8.70144992e+01, -1.06593964e+01,  1.29664573e+01,
                  -9.48925211e+00,  8.39629251e+01, -2.88142490e-01,  1.72892845e+00,
                   7.34808940e+01, -7.23146350e+00,  1.48900080e+00,  2.16863469e+00,
                  -5.94001900e-01,  9.59852500e-02,  1.09938000e-02, -1.90058600e-02,
                  -4.50952870e-01,  0.00000000e+00,  0.00000000e+00,  2.40009140e-01,
                   4.85970410e-01,  3.88285965e+00,  1.05845000e-01, -7.70040800e-02,
                   3.35122179e+00,  1.12624459e+01,  4.34711646e+00, -6.19191860e-01,
                   9.92973680e-01,  6.11405931e+00, -5.60836840e-01, -2.93954780e-01,
                   6.90423300e-02,  3.39794700e-02, -1.24912730e-01,  0.00000000e+00,
                   0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                  -4.85997170e-01,  3.89679000e-02, -1.31052600e-02,  5.76526800e-01,
                  -3.85070200e-02, -6.48371900e-01,  1.31382710e-01])

    # Issues with rounding
    assert np.allclose(coupling_array[0], coupling0, 0.00000001)

    couplinglen3 = [3, 2, 1, 0, 1, 2, 3, 4, 5, 0, 4, 0, 0, 5, 5, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 5, 4, 4, 5, 5, 2, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 0, 0]

    assert np.array_equal(coupling_len[3], couplinglen3)


def test_orcaget_opt_status():

    file = 'tests/test_store/TST_pass_orca_opt.log'
    status = orcaread.get_opt_status(file)
    assert status == 'successful'

    # Dont actually have a failed one yet
    '''
    file = 'tests/test_store/TST_fail_orca_opt.log'
    status = orcaread.get_opt_status(file)
    assert status == 'failed'
    '''

def test_orcaget_nmr_status():

    file = 'tests/test_store/TST_pass_orca_nmr.log'
    status = orcaread.get_nmr_status(file)
    assert status == 'successful'


    file = 'tests/test_store/TST_fail_orca_nmr.log'
    status = orcaread.get_nmr_status(file)
    assert status == 'failed'

def test_orcaread_opt():

    file = 'tests/test_store/TST_pass_orca_opt.log'
    energy = orcaread.read_opt(file)
    assert energy == -1072.556525704736

def test_orcaread_functional():

    file = 'tests/test_store/TST_pass_orca_nmr.log'
    functional, basisset = orcaread.read_functional(file)

    assert functional == 'wB97X-D3'
    assert basisset == '6-311g'

def test_orcaread_nmr():

    file = 'tests/test_store/TST_pass_orca_nmr.log'

    shift, coupling = orcaread.read_nmr(file)

    shiftVL = np.asarray([141.369, 123.42, 137.777, 157.912, 126.292, 167.17, 187.061, 175.828,  99.682,
                            269.125,  98.381, 177.817, 142.869,  40.027, 230.944,  38.064, 371.694, 183.907,
                            132.017, 220.953, 219.118, 206.043,  34.843,  31.793, 214.094, -52.769,  52.191,
                            110.627, 149.995,  37.441,  28.192,  32.913,  34.683,  29.322,  32.26,  31.577,
                            36.024,  27.327,  31.849,  29.532,  38.429,  34.514,  47.652,  34.18,  39.965,
                            32.893,  31.298])
    assert np.array_equal(shift, shiftVL)

    # Too many values to construct array so check a few
    assert coupling[0][0] == 0.0
    assert coupling[0][1] == 88.589
    assert coupling[3][0] == 13.592
    assert coupling[-1][-1] == 0.0
    assert coupling[6][-1] == -3.073
    assert coupling[4][34] == -0.01

def test_generic_pybel_read():
    files = ['tests/test_store/ethane.xyz', 'tests/test_store/ethane.mol', # 'tests/test_store/ethane_nmr.log',
            'tests/test_store/ethane.mol2', 'tests/test_store/ethane.sdf']
    intypes = ['xyz', 'mol',#'g09',
                'mol2', 'sdf']

    xyz, types, conn_table, coupling_len = strucread.generic_pybel_read(files[0], intypes[0])

    assert np.array_equal(xyz, [[-0.4125,  0.0000,  0.0000],
                     [ 0.4125,  0.0000,  0.0000],
                     [-0.9475,  0.9266,  0.0000],
                     [-0.9475, -0.9266,  0.0000],
                     [-1.4825, -0.0000,  0.0000],
                     [ 0.9475, -0.9266,  0.0000],
                     [ 0.9475,  0.9266,  0.0000],
                     [ 1.4825,  0.0000,  0.0000]])

    assert np.array_equal(types, [6, 6, 1, 1, 1, 1, 1, 1])

    assert np.array_equal(conn_table, [[0, 1, 1, 1, 1, 0, 0, 0],
                                         [1, 0, 0, 0, 0, 1, 1, 1],
                                         [1, 0, 0, 0, 0, 0, 0, 0],
                                         [1, 0, 0, 0, 0, 0, 0, 0],
                                         [1, 0, 0, 0, 0, 0, 0, 0],
                                         [0, 1, 0, 0, 0, 0, 0, 0],
                                         [0, 1, 0, 0, 0, 0, 0, 0],
                                         [0, 1, 0, 0, 0, 0, 0, 0]])

    assert np.array_equal(coupling_len, [[0, 1, 1, 1, 1, 2, 2, 2],
                                         [1, 0, 2, 2, 2, 1, 1, 1],
                                         [1, 2, 0, 2, 2, 3, 3, 3],
                                         [1, 2, 2, 0, 2, 3, 3, 3],
                                         [1, 2, 2, 2, 0, 3, 3, 3],
                                         [2, 1, 3, 3, 3, 0, 2, 2],
                                         [2, 1, 3, 3, 3, 2, 0, 2],
                                         [2, 1, 3, 3, 3, 2, 2, 0]])

    for file, type in zip(files, intypes):
        _xyz, _types, _conn_table, _coupling_len = strucread.generic_pybel_read(file, type)
        assert np.array_equal(xyz, _xyz)
        assert np.array_equal(types, _types)
        assert np.array_equal(conn_table, _conn_table)
        assert np.array_equal(coupling_len, _coupling_len)

def test_fast_generic_pybel_read():
    files = ['tests/test_store/ethane.xyz', 'tests/test_store/ethane.mol', # 'tests/test_store/ethane_nmr.log',
            'tests/test_store/ethane.mol2', 'tests/test_store/ethane.sdf']
    intypes = ['xyz', 'mol',#'g09',
                'mol2', 'sdf']

    xyz, types, conn_table, coupling_len = strucread.fast_generic_pybel_read(files[0], intypes[0])

    assert np.array_equal(xyz, [[-0.4125,  0.0000,  0.0000],
                     [ 0.4125,  0.0000,  0.0000],
                     [-0.9475,  0.9266,  0.0000],
                     [-0.9475, -0.9266,  0.0000],
                     [-1.4825, -0.0000,  0.0000],
                     [ 0.9475, -0.9266,  0.0000],
                     [ 0.9475,  0.9266,  0.0000],
                     [ 1.4825,  0.0000,  0.0000]])

    assert np.array_equal(types, [6, 6, 1, 1, 1, 1, 1, 1])

    assert np.array_equal(conn_table, [[0, 1, 1, 1, 1, 0, 0, 0],
                                         [1, 0, 0, 0, 0, 1, 1, 1],
                                         [1, 0, 0, 0, 0, 0, 0, 0],
                                         [1, 0, 0, 0, 0, 0, 0, 0],
                                         [1, 0, 0, 0, 0, 0, 0, 0],
                                         [0, 1, 0, 0, 0, 0, 0, 0],
                                         [0, 1, 0, 0, 0, 0, 0, 0],
                                         [0, 1, 0, 0, 0, 0, 0, 0]])

    assert np.array_equal(coupling_len, np.zeros((len(types), len(types)), dtype=np.float64))

    for file, type in zip(files, intypes):
        _xyz, _types, _conn_table, _coupling_len = strucread.fast_generic_pybel_read(file, type)
        assert np.array_equal(xyz, _xyz)
        assert np.array_equal(types, _types)
        assert np.array_equal(conn_table, _conn_table)
        assert np.array_equal(coupling_len, _coupling_len)


















 ##
