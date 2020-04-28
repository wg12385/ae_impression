import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')

from autoenrich.analysis.compare_mols import mol_isequal
from autoenrich.util.flag_handler.hdl_targetflag import target_to_flag
from autoenrich.molecule.nmrmol import nmrmol

from test_generators.dummy_mol import get_random_mol, get_ethane_mol

def test_mol_isequal():

    mol1 = get_ethane_mol()
    mol2 = get_random_mol()
    # Different molecules are different
    assert mol_isequal(mol1, mol2) == False
    # Identical molecules are identical
    assert mol_isequal(mol1, mol1) == True
