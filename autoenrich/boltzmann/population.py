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
import sys

# Get relative populations of conformers
def get_pop_array(conformers, temp=298):

    exclude = []

    for c, conformer in enumerate(conformers):
        if conformer.nmr_status != 'successful':
            exclude.append(0)
        else:
            exclude.append(1)

    if np.sum(exclude) == 0:
        print('No conformers to average, something went wrong. . .')
        sys.exit(0)

    e_array = np.zeros(len(conformers), dtype=np.float64)
    for c, conformer in enumerate(conformers):
        e_array[c] = conformer.energy

    kj_array = e_array * 2625.5
    min_val = np.amin(kj_array)
    rel_array = (kj_array - min_val) * exclude
    exp_array = -(rel_array*1000) / float(8.31*temp)
    exp_array = np.exp(exp_array) * exclude
    sum_val = np.sum(exp_array)

    pop_array = exp_array / sum_val

    return pop_array
