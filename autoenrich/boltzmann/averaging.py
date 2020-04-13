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

import numpy as np

# Botlzmann average the chemical shift values in a set of conformers
def boltzmann_shift(conformers):
	atoms = len(conformers[0].types)

	shift_array = np.zeros(atoms, dtype=np.float64)
	shift_var = np.zeros(atoms, dtype=np.float64)
	for conformer in conformers:
		if conformer.pop == 0:
			continue
		shift_array += conformer.shift * conformer.pop
		shift_var += conformer.shift_var * conformer.pop

	return shift_array, shift_var

# Boltzmann average the coupling array in a set of conformers
def boltzmann_coupling(conformers):
	atoms = len(conformers[0].types)

	coupling = np.zeros((atoms, atoms), dtype=np.float64)
	coupling_var = np.zeros((atoms, atoms), dtype=np.float64)
	for conformer in conformers:
		if conformer.pop == 0:
			continue
		coupling += conformer.coupling * conformer.pop
		coupling_var += conformer.coupling_var * conformer.pop

	return coupling, coupling_var
