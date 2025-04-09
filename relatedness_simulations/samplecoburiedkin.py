##Scripts for simulations of within-building co-burial relatedness in paper "Female lineages and changing kinship patterns in Neolithic Catalhoyuk"
##https://github.com/eren-yuncu

import pandas as pd
import numpy as np
from itertools import product
import random
import argparse

parser = argparse.ArgumentParser(description = 'Simulate co-burial relatedness')
parser.add_argument('-b', '--burial', type = str, required = True, help = 'Burial information from Table S1')
parser.add_argument('-k', '--kin', type = str, required = True, help = 'Observed kinship relatedness from Table S3')
parser.add_argument('-c', '--filter', type = int, required = True, help = 'Min common SNP cut off value for filtering observed kinship results')
parser.add_argument('-s', '--siter', type = int, required = True, help = 'Simulation no')
parser.add_argument('-t', '--type', type = str, required = True, help = 'Type of simulation: "cons" for constant size, "var" variable size')
parser.add_argument('-f', '--fam', type = int, help = 'Family size for constant size simulations')
parser.add_argument('-fsub', '--subfammin', type = int, required = True, help = 'Minimun family size for creating array if burial size is smaller than sample size')
parser.add_argument('-fmin', '--fammin', type = int, help = 'Min family size for variable size simulations')
parser.add_argument('-fmax', '--fammax', type = int, help = 'Max family size for variable size simulations')
args = parser.parse_args()

if args.type != 'cons' and args.type != 'var':
	print('Please input a valid simulation type: "cons" for constant size, "var" variable size')
	exit()

if args.type == 'cons':
	if args.fam is None:
		parser.error('Family size is required for constant family size simulations!')

if args.type == 'var':
	if args.fammin is None or args.fammax is None:
		parser.error('Min and max family size is required for variable family size simulations!')

if args.type == 'var':
	if args.fammin > args.fammax:
		parser.error('Min family size can not be bigger than max family size!')

#Read burial information
inds = pd.read_csv(args.burial, sep = '\t', converters = {'Building': str, 'Unit': str})

#Remove information of unmerged identical individuals
inds = inds[inds['Sample ID'] != 'cch328']
inds = inds[inds['Sample ID'] != 'cch320']
inds = inds[inds['Sample ID'] != 'cch341']
inds = inds[inds['Sample ID'] != 'cch348']

#Select and rename necessary columns
inds = inds.replace('N/A','').fillna('')

inds = inds[['Unit', 'Period', 'Building', 'Age', 'Kinship Analysis']]
inds.loc[(inds['Age'] == 'child/infant'), 'Age'] = 'infant'

#Assign age categories
sub = ['prenatal', 'neonate', 'infant', 'child']
adult = ['young adult', 'old adult', 'middle adult', 'adult']

inds.loc[inds.Age.isin(sub), 'Age_Cat'] = 'Subadult'
inds.loc[inds.Age.isin(adult), 'Age_Cat'] = 'Adult'
inds.loc[(inds.Age == 'adolescent'), 'Age_Cat'] = 'Adolescent'
inds.loc[(inds.Age == ''), 'Age_Cat'] = 'Unknown'
inds.loc[(inds.Age == 'age unknown'), 'Age_Cat'] = 'Unknown'

inds = inds[inds['Age_Cat'] != 'Unknown']

#Calculate number of burials in each building
inds['Period_Building'] = inds[['Period', 'Building']].agg('_'.join, axis=1)
inds['Period_BuildingjUnit'] = inds[['Period_Building', 'Unit']].agg('j'.join, axis=1)

unit_burials = pd.DataFrame({'Burials': inds.groupby(['Period_BuildingjUnit'])['Period_BuildingjUnit'].count()}).reset_index()
unit_burials[['Period_Building', 'Unit']] = unit_burials.Period_BuildingjUnit.str.split('j', expand = True)
burials = pd.DataFrame({'Burials': unit_burials.groupby(['Period_Building'])['Burials'].sum()}).reset_index()

age_table_sub = inds[['Period_BuildingjUnit', 'Age_Cat']]
age_table = age_table_sub.copy()
age_table[['Period_Building', 'Unit']] = age_table.Period_BuildingjUnit.str.split('j', expand = True)
age_table = age_table.drop(['Unit', 'Period_BuildingjUnit'], axis = 1)
age_freq = age_table.pivot_table(index = 'Period_Building', columns = 'Age_Cat', aggfunc = lambda x: len(x)).reset_index()

burials = pd.merge(burials, age_freq, on = 'Period_Building', how = 'outer')

#Calculate number of individuals with genetic info in each building
dataset = inds[inds['Kinship Analysis'] ==  'Yes']

age_table_temp = dataset[['Period_BuildingjUnit', 'Age_Cat']]
age_table_dna = age_table_temp.copy()
age_table_dna[['Period_Building', 'Unit']] = age_table_dna.Period_BuildingjUnit.str.split('j', expand = True)
age_table_dna = age_table_dna.drop(['Unit', 'Period_BuildingjUnit'], axis = 1)
age_freq_dna = age_table_dna.pivot_table(index = 'Period_Building', columns = 'Age_Cat', aggfunc = lambda x: len(x)).reset_index()

dna = pd.DataFrame({'DNA': dataset.groupby(['Period_Building'])['Period_Building'].count()}).reset_index()
dna = pd.merge(dna, age_freq_dna, on = 'Period_Building', how = 'outer')

buildings = list(dna['Period_Building'])
burials = burials[burials['Period_Building'].isin(buildings)].reset_index(drop = True)

#Merge burial and DNA information
data_all = pd.merge(burials, dna, on = 'Period_Building', how = 'outer')

data_all.columns = data_all.columns.str.replace('_x', '_Burials', regex = True)
data_all.columns = data_all.columns.str.replace('_y', '_DNA', regex = True)

data_all = data_all.fillna(0)
cols_list = data_all.columns.tolist()
cols_list.remove('Period_Building')

for i in cols_list:
	data_all[i] = data_all[i].astype(int)

data_all[['Period', 'Building']] = data_all.Period_Building.str.split('_', expand = True)
data_all = data_all.drop(['Period_Building'], axis = 1).iloc[:, [8,9,0,1,2,3,4,5,6,7]]

#Save burial and DNA information
data_all.to_csv('burial_DNA_info.txt', sep = '\t', index = False)

#Calculate number of coburial pairs with kinship info in each building
#Read kinsip results
res = pd.read_csv(args.kin, sep = '\t', converters = {'Building1': str, 'Building2': str}, low_memory = False)
res = res.replace('N/A','').fillna('')

#Rename necessary columns
res.columns = res.columns.str.replace('\n', '', regex=True)
res.columns = res.columns.str.replace(' ', '', regex=True)
res.columns = res.columns.str.replace('Ind', '', regex=True)

res.rename(columns={res.columns[3]: 'Relatedness_Auto'}, inplace = True)
res.rename(columns={res.columns[7]: 'CommonSNPs_READ2_Auto'}, inplace = True)
res.rename(columns={res.columns[13]: 'CommonSNPs_READ2_Imp_Auto'}, inplace = True)

#Select coburial pairs
res = res[res['Period1'] == res['Period2']]
res = res[res['Period1'] != '']
res = res[res['Building1'] == res['Building2']]
res = res[res['Building1'] != '']
res = res[res['Relatedness_Auto'] != '']

#Set SNP cut of value
res['CommonSNPs_READ2_Auto'] = pd.to_numeric(res['CommonSNPs_READ2_Auto'])
res['CommonSNPs_READ2_Imp_Auto'] = pd.to_numeric(res['CommonSNPs_READ2_Imp_Auto'])

if args.filter != 1000:	
	res = res[(res['CommonSNPs_READ2_Auto'] >= args.filter) | (res['CommonSNPs_READ2_Imp_Auto'] >= args.filter)]

#Select and rename necessary columns
res.loc[(res['Age1'] == 'child_infant'), 'Age1'] = 'infant'
res.loc[(res['Age2'] == 'child_infant'), 'Age2'] = 'infant'
res = res[(res['Period1'] == 'Early') | (res['Period1'] == 'Middle') | (res['Period1'] == 'Late')]

#Assign age categories
sub = ['prenatal', 'neonate', 'infant', 'child']
adult = ['young adult', 'old adult', 'middle adult', 'adult']

res.loc[res.Age1.isin(sub), 'Age1_Cat'] = 'Subadult'
res.loc[res.Age1.isin(adult), 'Age1_Cat'] = 'Adult'
res.loc[(res.Age1 == 'adolescent'), 'Age1_Cat'] = 'Adolescent'

res.loc[res.Age2.isin(sub), 'Age2_Cat'] = 'Subadult'
res.loc[res.Age2.isin(adult), 'Age2_Cat'] = 'Adult'
res.loc[(res.Age2 == 'adolescent'), 'Age2_Cat'] = 'Adolescent'

res['Period_Building'] = res[['Period1', 'Building1']].agg('_'.join, axis=1)
res['Age_Cat'] = ['_'.join(x) for x in np.sort(res[['Age1_Cat', 'Age2_Cat']], axis=1)]

#Calculate number of coburial pairs in each building
pairs = pd.crosstab(index=res['Period_Building'], columns=res['Age_Cat'])
pairs.reset_index(inplace=True)

pairs = pairs.fillna(0)
cols_list = pairs.columns.tolist()
cols_list.remove('Period_Building')

for i in cols_list:
	pairs[i] = pairs[i].astype(int)

pairs[['Period', 'Building']] = pairs.Period_Building.str.split('_', expand = True)
pairs = pairs.drop(['Period_Building'], axis = 1).iloc[:, [6,7,0,1,2,3,4,5]]

#Save coburial pair information
cut_off = str(args.filter)
pairs.to_csv('period_pair_info_filt' + cut_off + '.txt', sep = '\t', index = False)

#Start simulations
#Functions
def age_number(building):
	buildings = data[data['Period_Building'].isin(building)].reset_index(drop = True)
	age_table = buildings[['Adolescent_Burials', 'Adult_Burials', 'Subadult_Burials']]
	age_table = age_table.replace('', 0)
	age_table = age_table.astype(int)
	age_table.columns = age_table.columns.str.rstrip('_Burials') 
	age_no = age_table.to_dict('list')
	return age_no

#Filter data to use in simulations
data = data_all[(data_all['Period'] == 'Early') | (data_all['Period'] == 'Middle') | (data_all['Period'] == 'Late')]
data = data[data['Burials'] != 1]
data = data[data['DNA'] != 1]
data = data[data['DNA'] != 0]
data = data[data['Building'] != '']
data.dropna(subset = ['Building'], inplace = True)

data['Period_Building'] = data[['Period', 'Building']].agg('_'.join, axis=1)
buildings = list(data['Period_Building'])

#Assign individuals to buildings
temp1 = data.loc[data.index.repeat(data.Burials)].reset_index(drop = True)

col_names = list(temp1.columns)
col_names.append('Age')
temp = pd.DataFrame(columns = col_names)
for b in buildings:
	build = [b]
	groups = temp1[temp1['Period_Building'].isin(build)].reset_index(drop = True)
	temp_age = pd.DataFrame(columns = col_names)
	age_dict = age_number(build)
	for key, value in age_dict.items():
		temp_build = groups.sample(n = value[0])
		temp_build['Age'] = key
		temp_age = pd.concat([temp_age, temp_build])
	temp = pd.concat([temp, temp_age])

#Assign family size
if args.type == 'cons':
	temp = temp[['Period_Building', 'Burials', 'Age']]

	fam_no = args.fam
	sim = pd.DataFrame(columns = ['Period_Building', 'Burials', 'Age', 'Fam'])
	for b in buildings:
		build = [b]
		groups = temp[temp['Period_Building'].isin(build)].reset_index(drop = True)
		fams = []
		total = len(groups)
		f_list = groups.index.values.tolist()
		max_f = max([j for j in f_list if j % fam_no == 0])
		sub_fam_size = list(range(args.subfammin, (total + 1)))
		sub_fam_no = random.choice(sub_fam_size)
		sub_max_f = max([j for j in f_list if j % sub_fam_no == 0])
		for i, row in groups.iterrows():
			if total % fam_no == 0:
				fam = f'f{i // fam_no + 1}'
				fams.append(fam)
			else:
				if total < fam_no:
					if total % sub_fam_no == 0:
						fam = f'f{i // sub_fam_no + 1}'
					elif (i+1) > sub_max_f:
						fam = f'f{i}'
					else:
						fam = f'f{i // sub_fam_no + 1}'
				elif (i+1) > max_f:
					fam = f'f{i}'
				else:
					fam = f'f{i // fam_no + 1}'
				fams.append(fam)
		groups['Fam']  = fams
		sim = pd.concat([sim, groups])
	sim = sim.reset_index(drop = True)

if args.type == 'var':
	fam_size = list(range(args.fammin, (args.fammax + 1)))

	temp = temp[['Period_Building', 'Burials', 'Age']]

	sim = pd.DataFrame(columns = ['Period_Building', 'Burials', 'Age', 'Fam'])
	for b in buildings:
		build = [b]
		fam_no = random.choice(fam_size)

		groups = temp[temp['Period_Building'].isin(build)].reset_index(drop = True)
		fams = []
		total = len(groups)
		f_list = groups.index.values.tolist()
		max_f = max([j for j in f_list if j % fam_no == 0])
		sub_fam_size = list(range(args.subfammin, (total + 1)))
		sub_fam_no = random.choice(sub_fam_size)
		sub_max_f = max([j for j in f_list if j % sub_fam_no == 0])
		for i, row in groups.iterrows():
			if total % fam_no == 0:
				fam = f'f{i // fam_no + 1}'
				fams.append(fam)
			else:
				if total < fam_no:
					if total % sub_fam_no == 0:
						fam = f'f{i // sub_fam_no + 1}'
					elif (i+1) > sub_max_f:
						fam = f'f{i}'
					else:
						fam = f'f{i // sub_fam_no + 1}'
				elif (i+1) > max_f:
					fam = f'f{i}'
				else:
					fam = f'f{i // fam_no + 1}'
				fams.append(fam)
		groups['Fam']  = fams
		sim = pd.concat([sim, groups])
	sim = sim.reset_index(drop = True)

#Organize and save family table
sim['Family'] = sim['Period_Building'] + sim['Fam']
sim = sim.drop(['Burials', 'Fam'], axis = 1)

sim['Family'] = sim['Family'].str.replace('Early_', 'E')
sim['Family'] = sim['Family'].str.replace('Middle_', 'M')
sim['Family'] = sim['Family'].str.replace('Late_', 'L')

ind_no = pd.DataFrame({'ind_no': list(range(1, len(sim.index) +1))}).astype(str)
sim['Ind'] = sim['Family'] + ind_no['ind_no']
inds = sim['Ind'].to_list()

#Save family table
rln = sim.copy()
if args.type == 'cons':
	fam_s = str(args.fam)
else:
	fam_s = 'var'

siter_s = str(args.siter)
sim[['Period', 'Building']] = sim.Period_Building.str.split('_', expand = True)
sim = sim.drop(['Period_Building'], axis = 1).iloc[:, [3,4,2,1,0]]
sim.to_csv('simulation_main_distribution_fam' + fam_s + '_sim' + siter_s + '.txt', sep = '\t', index = False)

#Create relationship table
comb = pd.DataFrame(list(product(inds, inds)))
comb = comb[comb[0] != comb[1]]
comb[2] = ['-'.join(x) for x in np.sort(comb[[0, 1]], axis = 1)]
comb = comb.drop_duplicates(subset=[2], keep = 'first')
comb = comb.drop([2], axis = 1)

ind1 = pd.merge(comb, rln, left_on = 0, right_on = 'Ind', how = 'inner')
ind1 = ind1.rename(columns={0: 'Ind1'})
sim_rln = pd.merge(ind1, rln, left_on = 1, right_on = 'Ind', how = 'inner')
sim_rln = sim_rln.rename(columns={1: 'Ind2'})
sim_rln = sim_rln.drop(['Ind_x', 'Ind_y'], axis = 1)
sim_rln.columns = sim_rln.columns.str.replace('_x', '1', regex = True)
sim_rln.columns = sim_rln.columns.str.replace('_y', '2', regex = True)
sim_rln.loc[sim_rln['Family1'] == sim_rln['Family2'], 'Relatedness'] = 'Related'
sim_rln.loc[sim_rln['Family1'] != sim_rln['Family2'], 'Relatedness'] = 'Unrelated'

#Organize and save relationship table
sim_rln[['Period1', 'Building1']] = sim_rln.Period_Building1.str.split('_', expand = True)
sim_rln[['Period2', 'Building2']] = sim_rln.Period_Building2.str.split('_', expand = True)
sim_rln = sim_rln.drop(['Period_Building1', 'Period_Building2'], axis = 1).iloc[:, [7,8,0,3,2,9,10,1,5,4,6]]
sim_rln.to_csv('simulation_main_relatedness_fam' + fam_s + '_sim' + siter_s + '.txt', sep = '\t', index = False)

print('End of simulation!!')
