import json
import sys


def read_prefs(file):
	try:
		with open(file, 'r') as json_file:
			prefs = json.load(json_file)
	except Exception as E:
		print('Error reading preferences file ', file)
		print(E)
		print('Exiting. . .')
		sys.exit(0)


	return prefs

def write_default_prefs(file):

	prefs = {}
	prefs['conf_search'] = {}
	prefs['conf_search']['iterations'] = 200
	prefs['conf_search']['maxconfs'] = 100
	prefs['conf_search']['RMSthresh'] = 0.1
	prefs['conf_search']['Ethresh'] = 100000

	prefs['mol'] = {}
	prefs['mol']['charge'] = 0
	prefs['mol']['multiplicity'] = 1

	prefs['comp'] = {}
	prefs['comp']['parallel'] = True
	prefs['comp']['system'] = 'BC3'

	prefs['optimisation'] = {}
	prefs['optimisation']['memory'] = 26
	prefs['optimisation']['processors'] = 8
	prefs['optimisation']['opt'] = 'tight'
	prefs['optimisation']['freq'] = False
	prefs['optimisation']['functional'] = 'mpw1pw91'
	prefs['optimisation']['basisset'] = '6-311g(d,p)'
	prefs['optimisation']['solvent'] = 'none'
	prefs['optimisation']['solventmodel'] = 'iefpcm'
	prefs['optimisation']['grid'] = 'ultrafine'
	prefs['optimisation']['custom_cmd_line'] = False
	prefs['optimisation']['nodes'] = 1
	prefs['optimisation']['ppn'] = 8
	prefs['optimisation']['mem'] = 32
	prefs['optimisation']['walltime'] = '100:00:00'

	prefs['NMR'] = {}
	prefs['NMR']['memory'] = 26
	prefs['NMR']['processors'] = 8
	prefs['NMR']['mixed'] = True
	prefs['NMR']['functional'] = 'wb97xd'
	prefs['NMR']['basisset'] = '6-311g(d,p)'
	prefs['NMR']['solvent'] = 'none'
	prefs['NMR']['solventmodel'] = 'iefpcm'
	prefs['NMR']['custom_cmd_line'] = False
	prefs['NMR']['nodes'] = 1
	prefs['NMR']['ppn'] = 8
	prefs['NMR']['mem'] = 32
	prefs['NMR']['walltime'] = '100:00:00'

	json.dump(prefs, open(file, 'w'), indent=4)
