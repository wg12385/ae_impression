# Copyright 2020 Will Gerrard
#This file is part of autoENRICH.

#autoENRICH is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoENRICH is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with autoENRICH.  If not, see <https://www.gnu.org/licenses/>.

import pybel as pyb
from autoENRICH.pybel_helpers import pybel_bonds, pybel_read
import numpy as np


def generic_pybel_read(file, type):
	mol = next(pyb.readfile(type, file))
	type_list, types = pybel_read.mol_read_type(mol)
	xyz = pybel_read.mol_read_xyz(mol)

	conn_table = pybel_bonds.mol_get_bond_table(mol)

	coupling_len = pybel_bonds.get_coupling_lengths(mol, types, maxlen=6)

	return xyz, types, conn_table, coupling_len


def fast_generic_pybel_read(file, type):
	mol = next(pyb.readfile(type, file))
	type_list, types = pybel_read.mol_read_type(mol)
	xyz = pybel_read.mol_read_xyz(mol)

	conn_table = pybel_bonds.mol_get_bond_table(mol)

	#coupling_len = pybel_bonds.get_coupling_lengths(mol, types, maxlen=6)
	coupling_len = np.zeros((len(types), len(types)), dtype=np.float64)

	return xyz, types, conn_table, coupling_len
