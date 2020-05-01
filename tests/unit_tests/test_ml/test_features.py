# Copyright 2020 Will Gerrard
#This file is part of autoenrich.

#autoenrich is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoenrich is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with autoenrich.  If not, see <https://www.gnu.org/licenses/>.
import sys
import os

import numpy as np

from autoenrich.ml.features.GNR_features import get_dummy_features

from test_generators.dummy_mol import get_ethane_mol


def test_get_dummy_features():

    mol = get_ethane_mol()

    mols = [mol]

    x, y, r = get_dummy_features(mols, targetflag='CCS')
    assert x == []
    assert np.array_equal(y, mol.shift[np.where(mol.types==6)])
    assert np.array_equal(r, [['ethane', 0], ['ethane', 1]])

    x, y, r = get_dummy_features(mols, targetflag='HCS')
    assert x == []
    assert np.array_equal(y, mol.shift[np.where(mol.types==1)])
    assert np.array_equal(r, [['ethane', 2], ['ethane', 3], ['ethane', 4],
                            ['ethane', 5], ['ethane', 6], ['ethane', 7]])

    x, y, r = get_dummy_features(mols, targetflag='1JCH')
    for i, index in enumerate(r):
        assert y[i] == mol.coupling[index[1]][index[2]]
    assert np.array_equal(r, [['ethane', 0, 2], ['ethane', 0, 3], ['ethane', 0, 4],
                            ['ethane', 1, 5], ['ethane', 1, 6], ['ethane', 1, 7]])

    x, y, r = get_dummy_features(mols, targetflag='3JHH')
    for i, index in enumerate(r):
        assert y[i] == mol.coupling[index[1]][index[2]]
    assert np.array_equal(r, [['ethane', 2, 5], ['ethane', 2, 6], ['ethane', 2, 7],
                            ['ethane', 3, 5], ['ethane', 3, 6], ['ethane', 3, 7],
                            ['ethane', 4, 5], ['ethane', 4, 6], ['ethane', 4, 7]])


def test_qml_features():

    import qml

    mol = get_ethane_mol()
    mbtypes = [[1],[1,1], [1,1,1], [1,1,6], [1,1,7], [1,1,8], [1,1,9], [1,6], [1,6,1], [1,6,6], [1,6,7], [1,6,8], [1,6,9], [1,7], [1,7,1], [1,7,6], [1,7,7], [1,7,8], [1,7,9], [1,8], [1,8,1], [1,8,6], [1,8,7], [1,8,8], [1,8,9], [1,9], [1,9,1], [1,9,6], [1,9,7], [1,9,8], [1,9,9], [6], [6,1], [6,1,1], [6,1,6], [6,1,7], [6,1,8], [6,1,9], [6,6], [6,6,1], [6,6,6], [6,6,7], [6,6,8], [6,6,9], [6,7], [6,7,1], [6,7,6], [6,7,7], [6,7,8], [6,7,9], [6,8], [6,8,1], [6,8,6], [6,8,7], [6,8,8], [6,8,9], [6,9], [6,9,1], [6,9,6], [6,9,7], [6,9,8], [6,9,9], [7],[7,1], [7,1,1], [7,1,6], [7,1,7], [7,1,8], [7,1,9], [7,6], [7,6,1], [7,6,6], [7,6,7], [7,6,8], [7,6,9], [7,7], [7,7,1], [7,7,6], [7,7,7], [7,7,8], [7,7,9], [7,8], [7,8,1], [7,8,6], [7,8,7], [7,8,8], [7,8,9], [7,9], [7,9,1], [7,9,6], [7,9,7], [7,9,8], [7,9,9], [8], [8,1], [8,1,1], [8,1,6], [8,1,7], [8,1,8], [8,1,9], [8,6], [8,6,1], [8,6,6], [8,6,7], [8,6,8], [8,6,9], [8,7], [8,7,1], [8,7,6], [8,7,7], [8,7,8], [8,7,9], [8,8], [8,8,1], [8,8,6], [8,8,7], [8,8,8], [8,8,9], [8,9], [8,9,1], [8,9,6], [8,9,7], [8,9,8], [8,9,9], [9], [9,1], [9,1,1], [9,1,6], [9,1,7], [9,1,8], [9,1,9], [9,6], [9,6,1], [9,6,6], [9,6,7], [9,6,8], [9,6,9], [9,7], [9,7,1], [9,7,6], [9,7,7], [9,7,8], [9,7,9], [9,8], [9,8,1], [9,8,6], [9,8,7], [9,8,8], [9,8,9], [9,9], [9,9,1], [9,9,6], [9,9,7], [9,9,8], [9,9,9]]
    reps = qml.representations.generate_slatm(mol.xyz, mol.types, mbtypes)
    x = np.asarray(reps)
    reps = qml.representations.generate_atomic_coulomb_matrix(mol.types, mol.xyz)
    x = np.asarray(reps)
    reps = qml.fchl.generate_representation(mol.xyz, mol.types)
    x = np.asarray(reps)
    try:
        reps = qml.representations.generate_acsf(mol.types, mol.xyz)
    except:
        print('Warning: not using qml version with ACSF rep')











###
