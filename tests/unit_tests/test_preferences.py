


import os

from autoenrich.preferences.preferences import read_prefs, check_prefs, get_default_prefs, write_default_prefs


def test_read_prefs():

    file = 'tests/test_store/ENRICH.json'

    prefs = read_prefs(file)

    assert isinstance(prefs, dict)
    for pref in prefs.values():
        assert isinstance(prefs, dict)

def test_get_default_prefs():

    prefs = get_default_prefs()
    assert isinstance(prefs, dict)
    for pref in prefs.values():
        assert isinstance(prefs, dict)

def test_check_prefs():

    prefs = get_default_prefs()

    prefs, chg = check_prefs(prefs)
    assert not chg

    del prefs['NMR']
    del prefs['optimisation']['functional']

    prefs, chg = check_prefs(prefs)
    assert chg
    assert 'NMR' in prefs.keys()
    assert 'functional' in prefs['optimisation'].keys()


def test_write_default_prefs():

    file = 'tests/test_tmp/prefs.tmp'

    write_default_prefs(file)

    assert os.path.isfile(file)

    prefs = read_prefs(file)
    assert isinstance(prefs, dict)
    for pref in prefs.values():
        assert isinstance(prefs, dict)

    os.remove(file)















#
