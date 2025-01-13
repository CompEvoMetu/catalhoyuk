##Scripts for simulations of within-building co-burial relatedness in paper "Female lineages and changing kinship patterns in Neolithic Catalhoyuk"
##https://github.com/eren-yuncu

import pandas as pd
import numpy as np
import glob, os
import warnings
from scipy.stats import spearmanr
import argparse

parser = argparse.ArgumentParser(description = 'Calculate simulated rho')
parser.add_argument('-f', '--fam', type = int, help = 'Family size for constant size simulations')
parser.add_argument('-t', '--type', type = str, required = True, help = 'Type of simulation: "cons" for constant size, "var" variable size')
parser.add_argument('-i', '--iter', type = int, required = True, help = 'Total number of iterations for randomization')
parser.add_argument('-c', '--filter', type = int, required = True, help = 'Min common SNP cut of value for filtering observed kinship results')
args = parser.parse_args()

if args.type != 'cons' and args.type != 'var':
	print('Please input a valid simulation type: "cons" for constant size, "var" variable size')
	exit()

if args.type == 'cons':
	if args.fam is None:
		parser.error('Family size is required for constant family size simulations!')

if args.type == 'cons':
	fam_s = str(args.fam)
else:
	fam_s = 'var'

cut_of = str(args.filter)
iter_s = str(args.iter)

warnings.filterwarnings('ignore', message = 'An input array is constant; the correlation coefficient is not defined.')

#Concat results of each simulation
results = glob.glob('simulation_random_relatedness_per_building_freq*fam' + fam_s + '_sim*_filt' + cut_of + '.txt')

if len(results) == 0:
	print('Can not find simulated frequency values for buildings. Please run samplecobureidkin.py and randomizecoburiedkin.py first or check parameters!')
	exit()

#Get column names
temp = pd.read_csv(results[0], sep = '\t')
col_names = list(temp.columns.values)

#List ages
sub1 = '_relatedness'
sub2 = '_fam'

ages = []
for i in results:
	idx1 = i.index(sub1)
	idx2 = i.index(sub2)
	x = ''
	for idx in range(idx1 + len(sub1) + 1, idx2):
		x = x + i[idx]
		age = x + '_fam'
	ages.append(age)

age_list = list(np.unique(ages))

for a in age_list:
	merge_list = []
	for result in results:
		if a in result:
			merge_list.append(result)

	#Merge different simulations and calculate rho
	df = pd.DataFrame(columns = col_names)
	for name in merge_list:
		res = pd.read_csv(name, sep = '\t')
		df = pd.concat([df, res])
		df['Simulation'] = df['Simulation'].astype(int)
		df['Iteration'] = df['Iteration'].astype(int)
		df['Unrelated'] = df['Unrelated'].astype(int)
		df['Related'] = df['Related'].astype(int)
		df['Total'] = df['Total'].astype(int)
		df['Simulation'] = df['Simulation'].astype(str)
		df['Iteration'] = df['Iteration'].astype(str)
		#Save merged results
		df.to_csv('simulation_random_relatedness_' + a + fam_s + '_iter'+ iter_s + '_filt' + cut_of + '.txt', sep = '\t', index = False)

		#Get building date information
		df[['Min', 'Max']] = df.Date.str.split('-', expand = True)
		df['Min'] = df['Min'].astype(int)
		df['Max'] = df['Max'].astype(int)
		df['Mean'] = (df['Min'] + df['Max'])/2

		#Filter buildings with < 2 burials
		df['Prop'] = df['Related']/df['Total']
		df = df[df['Total'] > 2]

		df['Groups'] = df[['Simulation', 'Iteration']].agg('_'.join, axis=1)

		iters = df['Groups'].unique().tolist()

		rho_all = pd.DataFrame(columns = ['rho', 'p_Value', 'Groups'])

		for i in iters:
			df_a = df[df['Groups'] == i].reset_index(drop = True)
			rho, p_val = spearmanr(df_a['Prop'], df_a['Mean'])

			rho_res = pd.DataFrame(columns = ['rho', 'p_Value', 'Groups'])
			rho_res['rho'] = [rho]
			rho_res['p_Value'] = [p_val]
			rho_res['Groups'] = [i]
			rho_all = pd.concat([rho_all, rho_res])

		rho_all[['Simulation', 'Iteration']] = rho_all['Groups'].str.split('_', expand = True)
		rho_all = rho_all.drop('Groups', axis = 1).iloc[:, [2,3,0,1]]

		rho_all.to_csv('rho_burial2_simulation_random_relatedness_' + fam_s + '_iter' + iter_s + '_filt' + cut_of + '.txt', na_rep = 'NaN', sep = '\t', index = False)

print('End of rho calculation!!')
