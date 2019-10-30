import pybel as pyb
from pybel_helpers import pybel_bonds, pybel_read


def generic_pybel_read(file, type):
	mol = next(pyb.readfile(type, file))
	type_list, types = pybel_read.mol_read_type(mol)
	xyz = pybel_read.mol_read_xyz(mol)

	conn_table = pybel_bonds.mol_get_bond_table(mol)

	coupling_len = pybel_bonds.get_coupling_lengths(mol, types, maxlen=6)

	return xyz, types, dist, conn_table, coupling_len
