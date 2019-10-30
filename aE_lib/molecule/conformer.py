import numpy as np
import glob
from .nmrmol import nmrmol

class conformer(nmrmol):

	def __init__(self, molid, path=''):
		nmrmol.__init__(self, path=path, molid=molid)

		# store location of conf search xyz file
		self.xyz_file = 'None'

		# store location and status of optimisation files
		self.opt_in = 'None'
		self.opt_log = 'None'
		self.opt_status = 'None'

		# store location and status of NMR files
		self.nmr_in = 'None'
		self.nmr_log = 'None'
		self.nmr_status = 'None'

		self.energy = 404.404
		self.pop = 404.404




###
