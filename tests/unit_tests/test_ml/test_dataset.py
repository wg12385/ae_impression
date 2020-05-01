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
import glob

from autoenrich.molecule.dataset import dataset
from autoenrich.file_creation.structure_formats.nmredata import nmrmol_to_nmredata

from test_generators.dummy_mol import get_random_ethane, get_random_mol, mol_is_ethane



def test_dataset():

    dset = dataset()

    assert dset.mols == []
    assert dset.x == []
    assert dset.y == []
    assert dset.r == []
    assert dset.mol_order == []
    assert dset.big_data == False
    assert dset.files == []
    assert dset.params == {}

def test_get_mols():

    files = glob.glob('tests/test_store/dataset/eth_rnd_mol_*.nmredata.sdf')
    dset = dataset()
    dset.get_mols(files)

    assert dset.files == files
    assert dset.big_data == False
    assert len(dset.mols) == 5
    for mol in dset.mols:
        assert mol_is_ethane(mol)

'''
def test_get_features_from_mols():

    dset = dataset()
    for _ in range(5):
        dset.mols.append(get_random_ethane())



    dset.get_features_frommols(args)
'''















#
