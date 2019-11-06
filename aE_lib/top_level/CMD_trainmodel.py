
from IMPRESSION.models import FCHLmodel, KRRmodel, TFMmodel
import bayes_opt as bayes
from bayes_opt import BayesianOptimization
from bayes_opt import UtilityFunction
from bayes_opt.observer import JSONLogger
from bayes_opt.event import Events
from bayes_opt.util import load_logs

def train_model(x, y, modelflag='FCHL', featureflag='FCHL', searchflag='grid',
 				id='test_model', param_ranges={}):

	if searchflag == 'grid':
		score = grid_search(x, y, modelflag, featureflag, id, param_ranges)

	elif searchflag == 'gaussian':
		score = gaussian_search(x, y, modelflag, featureflag, id, param_ranges)


	print('Model optimised, score = ', score)
