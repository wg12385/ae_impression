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
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')

import numpy as np

from autoenrich.ml.features import QML_features

from test_generators.dummy_mol import get_ethane_mol


def test_get_FCHL_features():

    mol = get_ethane_mol()

    mols = [mol]

    x, y, r = get_dummy_features(mols, targetflag='CCS')
    assert x == []
    assert np.array_equal(y, mol.shift[np.where(mol.types==6)])
