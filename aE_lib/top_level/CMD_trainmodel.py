
from file_creation.HPS_submission import make_HPC_header
from IMPRESSION.train_model import HPS
from molecule.dataset import dataset
from molecule.nmrmol import nmrmol
import sys

# creates a submission file that will run a hyper parameter search
def setup_trainmodel(args):

	# NEED TO MAKE GPU VERSION
	python_env = args.python_env
	system = args.system
	nodes = 1
	ppn = args.processors
	walltime = args.walltime
	mem = args.memory

	jobname = 'IMP_' + args.modelflag + '_' + args.featureflag + '_' + args.targetflag + '_' + args.searchmethod + '_HPS'
	header = make_HPC_header(jobname, system=args.system, nodes=1, ppn=args.processors, walltime=args.walltime, mem=args.memory)
	strings.append('')
	strings.append('source activate {env:<10s}'.format(env=python_env))
	strings.append('IMPRESSION train --modelflag {model:<10s}  --featureflag {feature:<10s} \ '.format(model=args.modelflag,
																									feature=args.featureflag))
	strings.append('				 --paramflag {param:<10s}  --searchmethod {search:<10s} \ '.format(param=args.targetflag,
																									search=args.searchmethod))
	strings.append('				 --paramranges {param:<10s} \ '.format(param=args.param_ranges))
	strings.append('				 --paramlogs {param:<10s} \ '.format(param=args.param_logs))
	strings.append('				 --cv_steps {cv_steps:<10s}\ '.format(cv_steps=args.cv_steps))

	if args.searchmethod == 'grid':
		strings.append('				 --grid_density {grid_density:<10s}\ '.format(grid_density=args.grid_density))

	elif args.searchmethod == 'gaussian':
		strings.append('				 --epochs {epochs:<10s}\ '.format(epochs=args.epochs))
		strings.append('				 --kappa {kappa:<10s}\ '.format(kappa=args.kappa))
		strings.append('				 --xi {xi:<10s}\ '.format(xi=args.xi))
		strings.append('				 --random {random:<10s}\ '.format(random=args.random))

	elif args.searchmethod == 'random':
		strings.append('				 --epochs {epochs:<10s}\ '.format(epochs=args.epochs))


	with open(submission_file, 'w') as f:
		for string in header:
			print(string, file=f)
		for string in strings:
			print(string, file=f)


def train_model(args):

	print(args.training_set)
	sys.exit(0)


	if len(args.training_set) == 1:
		if args.training_set.split('.')[-1] == 'pkl':
			dset = pickle.load(open(args.training_set, 'rb'))
		elif args.training_set.split('.')[-1] == 'csv':
			dset = load_dataset_from_csv(args.training_set)
	else:
		files = glob.glob(args.training_set)
		dset = dataset()
		if args.store_datasets:
			dataset.get_mols(files, type='nmredata')
			dataset.get_features_frommols(files, args.featureflag, args.targetflag, args.cutoff)
			pickle.dump(dataset, open('training_set.pkl', 'wb'))
		else:
			dataset.get_features_fromfiles(files, args.featureflag, args.targetflag, args.cutoff)


	HPS(dataset.x, dataset.y, modelflag=args.modelflag, featureflag=args.featureflag, searchflag=args.searchmethod,
 				id='test_model', logfile=args.logfile, param_ranges=args.param_ranges,
				param_logs=args.param_logs, grid_density=args.grid_density, cv_steps=args.cv_steps,
				epochs=args.epochs, kappa=args.kappa, xi=args.xi, random=args.random)





###
