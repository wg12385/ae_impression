import numpy as np
import glob
from conversion.convert import *
from conversion.read_orca import eread, check_opt_status, check_nmr_status
from conversion.make_orca import make_optin, make_nmrin, make_submission_array, make_submission_qsub

class conformer(nmrmol):

	def __init__(self):
		nmrmol.__init(self))

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
