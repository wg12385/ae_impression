
from IMPRESSION.models.model import genericmodel

class FCHLmodel(genericmodel):

	def __init__(self, id='FCHLmodel', x=[], y=[], params={}):
		genericmodel.__init__(self, id, x, y, params)

		for key in ['sig', 'lam', 'cut']
