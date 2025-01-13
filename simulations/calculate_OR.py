##Scripts for simulations of within-building co-burial relatedness in paper "Female lineages and changing kinship patterns in Neolithic Catalhoyuk"
##https://github.com/eren-yuncu

import pandas as pd
import numpy as np
import glob, os
from scipy.stats import fisher_exact
import math
import argparse

parser = argparse.ArgumentParser(description = 'Calculate simulated OR')
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

#Concat results of each simulation
results = glob.glob('simulation_random_relatedness_freq*fam' + fam_s + '_sim*_filt' + cut_of + '.txt')

if len(results) == 0:
	print('Can not find simulated frequency values for periods. Please run samplecobureidkin.py and randomizecoburiedkin.py first or check parameters!')
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

#Merge different simulations and calculate OR
	df = pd.DataFrame(columns = col_names)
	for name in merge_list:
		res = pd.read_csv(name, sep = '\t', converters = {'Building': str}, low_memory = False)
		df = pd.concat([df, res])
		df['Simulation'] = df['Simulation'].astype(int)
		df['Iteration'] = df['Iteration'].astype(int)
		df['Period'] = pd.Categorical(df.Period, ordered = True, categories = ['Early', 'Middle', 'Late'])
		df = df.sort_values('Period')
		df = df.sort_values(['Simulation', 'Iteration'], ascending=[True, True])
		df['Unrelated'] = df['Unrelated'].astype(int).fillna(value = 0)
		df['Related'] = df['Related'].astype(int).fillna(value = 0)
		df['Total'] = df['Total'].astype(int).fillna(value = 0)
		df['Simulation'] = df['Simulation'].astype(str)
		df['Iteration'] = df['Iteration'].astype(str)
		#Save merged results
		df.to_csv('simulation_random_relatedness_' + a + fam_s + '_iter' + iter_s + '_filt' + cut_of + '.txt', sep = '\t', index = False)

		#Calculate OR
		
		df['Groups'] = df[['Simulation', 'Iteration']].agg('_'.join, axis=1)

		df['Related'] = df['Related'] + 1
		df['Unrelated'] = df['Unrelated'] + 1

		iters = df['Groups'].unique().tolist()

		or_all = pd.DataFrame(columns = ['Period', 'OR', 'logOR', 'Groups'])

		for i in iters:
			df_a = df[df['Groups'] == i].reset_index(drop = True)

			df_em = df_a[(df_a['Period'] == 'Early') | (df_a['Period'] == 'Middle')].reset_index(drop = True)
			df_em['Period'] = pd.Categorical(df_em.Period, ordered = True, categories = ['Early', 'Middle'])
			df_em = df_em.sort_values('Period')
			df_em = df_em[['Related', 'Unrelated']]
			or_em, p_em = fisher_exact(df_em, alternative = 'greater')
			lor_em = math.log(or_em)

			df_ml = df_a[(df_a['Period'] == 'Middle') | (df_a['Period'] == 'Late')].reset_index(drop = True)
			df_ml['Period'] = pd.Categorical(df_ml.Period, ordered = True, categories = ['Middle', 'Late'])
			df_ml = df_ml.sort_values('Period')
			df_ml = df_ml[['Related', 'Unrelated']]
			or_ml, p_ml = fisher_exact(df_ml, alternative = 'greater')
			lor_ml = math.log(or_ml)

			df_el = df_a[(df_a['Period'] == 'Early') | (df_a['Period'] == 'Late')].reset_index(drop = True)
			df_el['Period'] = pd.Categorical(df_el.Period, ordered = True, categories = ['Early', 'Late'])
			df_el = df_el.sort_values('Period')
			df_el = df_el[['Related', 'Unrelated']]
			or_el, p_el = fisher_exact(df_el, alternative = 'greater')
			lor_el = math.log(or_el)

			or_res = pd.DataFrame(columns = ['Period', 'OR', 'logOR', 'Groups'])
			or_res['Period'] = ['Early-Middle', 'Middle-Late', 'Early-Late']
			or_res['OR'] = [or_em, or_ml, or_el]
			or_res['logOR'] = [lor_em, lor_ml, lor_el]
			or_res['Groups'] = [i, i, i]
			or_all = pd.concat([or_all, or_res])

		or_all[['Simulation', 'Iteration']] = or_all['Groups'].str.split('_', expand = True)
		or_all['Simulation'] = or_all['Simulation'].astype(int)
		or_all['Iteration'] = or_all['Iteration'].astype(int)
		or_all = or_all.drop('Groups', axis = 1).iloc[:, [3,4,0,1,2]]

		#Save OR results
		or_all.to_csv('OR_simulation_random_relatedness_' + a + fam_s + '_iter' + iter_s + '_filt' + cut_of + '.txt', sep = '\t', index = False)

print('End of OR calculation!!')
