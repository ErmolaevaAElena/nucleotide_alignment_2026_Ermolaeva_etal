from modeller import * 

from modeller.automodel import * # Load the AutoModel class 

import sys 

log.verbose() 

# Override the 'special_restraints' and 'user_after_single_model' methods: 

class MyModel(AutoModel): 

	def user_after_single_model(self): 

# Report on symmetry violations greater than 1A after building

		self.restraints.symmetry.report(1.0)

	def special_patches(self, aln):

		self.rename_segments(segment_ids=['A'],renumber_residues=[1])

for pdbid in ['SCDN']:	

	env = Environ() 

# directories for input atom files

	env.io.atom_files_directory = ['.', '']

	a = MyModel(env,

		alnfile = '', # alignment filename

		knowns = ['7C2L', '7L2C', '7L2D', '7UKL', '8HLD'], # codes of the templates

		sequence = pdbid, # code of the target

		assess_methods=(assess.DOPE, assess.GA341)) 

	a.starting_model= 1

	a.ending_model = 3

	a.make()

# Get a list of all successfully built models from a.outputs

	ok_models = [x for x in a.outputs if x['failure'] is None]

# Rank the models by DOPE score

	key = 'DOPE score'

	if sys.version_info[:2] == (2,3):

		# Python 2.3's sort doesn't have a 'key' argument 

		ok_models.sort(lambda a,b: cmp(a[key], b[key]))

	else:

		ok_models.sort(key=lambda a: a[key])

# Get top model

	m = ok_models[0]

	print('Top model: %s (DOPE score %.3f)' % (m['name'], m[key]))

	ouf = open('%s_scores.csv' % pdbid, 'w')

	ouf.write('model_name,DOPE_score\n')

	for _model in ok_models:

		ouf.write('%s,%.3f\n' % (_model['name'], _model[key]))

	ouf.close()
