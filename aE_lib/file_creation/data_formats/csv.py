# Copyright 2020 Will Gerrard
#This file is part of autoENRICH.

#autoENRICH is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoENRICH is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with autoENRICH.  If not, see <https://www.gnu.org/licenses/>.

# write functions to output data in csv format ??
from reference.periodic_table import Get_periodic_table



def print_mol_csv(outname, refs, typerefs, values, labels):

	lines = []
	p_table = Get_periodic_table()

	sets = len(refs[0])

	idstring = ""
	for x in range(len(refs[0])):
		idstring = idstring + "{label:<s}MOLID,".format(label=labels[x])

	refstring = ""
	for y in range(1, len(refs[0][0])):
		refstring = refstring + "{atom:<s},{type:<s},".format(atom='Atom',
											type='Type',)

	valstring = ""
	for z in range(len(values[0])):
		valstring = valstring + "{label:<s}VALUE,".format(label=labels[z])

	lines.append(idstring+refstring+valstring)

	for i in range(len(refs)):


		idstring = ""
		for x in range(len(refs[i])):
			idstring = idstring + "{id:<s},".format(id=refs[i][x][0])

		refstring = ""
		for y in range(1, len(refs[0][0])):
			refstring = refstring + "{atom:<s},{type:<s},".format(atom=str(refs[i][0][y]),
												type=p_table[typerefs[i][0][y-1]])

		valstring = ""
		for z in range(len(values[i])):
			valstring = valstring + "{value:<.4f},".format(value=values[i][z])

		lines.append(idstring+refstring+valstring)


	with open(outname, 'w') as f:
		for line in lines:
			print(line, file=f)
