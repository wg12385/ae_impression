import sys
import os
sys.path.append(os.path.realpath(os.path.dirname(__file__))+'/../')


from autoenrich.file_creation.structure_formats.nmredata import nmrmol_to_nmredata
from autoenrich.file_creation.structure_formats.xyz import nmrmol_to_xyz
from autoenrich.file_creation.data_formats.csv import print_mol_csv
from autoenrich.file_creation.confsearch import make_confsearch_script
from autoenrich.file_creation.HPC_submission import get_chunks, make_HPC_header, make_HPC_batch_submission

from autoenrich.file_creation import orca_submission as orcasub
from autoenrich.file_creation import g09_submission as g09sub
from autoenrich.file_creation import g16_submission as g16sub

from test_generators.dummy_dataset import get_dummy_dataset
from test_generators.dummy_prefs import get_dummy_prefs
from test_generators.dummy_mol import get_random_mol


def test_csv():

    # Not written yet, weird test case, function not used much
    print('Not written yet, weird test case, function not used much')

def test_nmredata():

    mol = get_random_mol()

    outfile = 'tmp.nmredata.sdf'

    print(mol.coupling_len)

    nmrmol_to_nmredata(mol, outfile)
    # Check file was created
    assert os.path.exists(outfile)
    os.remove(outfile)

def test_xyz():

    mol = get_random_mol()

    outfile = 'tmp.xyz'

    nmrmol_to_xyz(mol, outfile)
    # Check file was created
    assert os.path.exists(outfile)
    os.remove(outfile)

def test_confsearch():

    scriptname = 'tmp.py'
    make_confsearch_script(scriptname, 'some_pickle.pkl', 'CCO', path='', iterations=2000,
    							RMSthresh=10.0, maxconfs=10, Ethresh=999.99)
    # Check file was created
    assert os.path.exists(scriptname)
    os.remove(scriptname)

def test_get_chunks():

    chunks = get_chunks(203)
    assert chunks == 5

def test_make_HPC_header():

    strings = make_HPC_header(jobname='auto-ENRICH',
                        system='PBS', nodes=1, ppn=1,
                            walltime="100:00:00", mem=3)

    # Nothing to test really . . .

def test_make_HPC_batch_submission():

    prefs = get_dummy_prefs()
    molname = 'dummy_mol'
    in_array = 'IN_ARRAY.txt'
    softwares = ['orca', 'g09', 'g16']
    start = 0
    end = 2
    for software in softwares:
        strings = make_HPC_batch_submission(prefs, molname, in_array, start, end,
                                                    software=software ,jobname='auto-ENRICH',
                                                        nodes=1, ppn=1, mem=3, walltime="100:00:00")

    # Nothing to test really

def test_make_optin():

    prefs = get_dummy_prefs()
    mol = get_random_mol()
    molname = 'dummy_mol'

    infile = orcasub.make_optin(prefs, molname, mol.xyz, mol.types)
    # Check file was created
    assert os.path.isfile(infile)
    os.remove(infile)

    infile = g09sub.make_optcom(prefs, molname, mol.xyz, mol.types)
    # Check file was created
    #assert os.path.isfile(infile)
    #os.remove(infile)
    infile = g16sub.make_optcom(prefs, molname, mol.xyz, mol.types)
    # Check file was created
    #assert os.path.isfile(infile)
    #os.remove(infile)














####
