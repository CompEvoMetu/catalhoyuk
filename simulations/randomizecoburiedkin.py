##Scripts for simulations of within-building co-burial relatedness in paper "Female lineages and changing kinship patterns in Neolithic Catalhoyuk"
##https://github.com/eren-yuncu

import numpy as np
import pandas as pd
from itertools import product
import argparse
import glob, os

parser = argparse.ArgumentParser(description = 'Randomize co-burial relatedness')
parser.add_argument('-d', '--date', type = str, required = True, help = 'Building date information from Table S26')
parser.add_argument('-i', '--iter', type = int, required = True, help = 'Number of iterations for randomization')
parser.add_argument('-s', '--siter', type = int, required = True, help = 'Simulation no')
parser.add_argument('-t', '--type', type = str, required = True, help = 'Type of simulation: "cons" for constant size, "var" variable size')
parser.add_argument('-f', '--fam', type = int, help = 'Family size for constant size simulations')
parser.add_argument('-c', '--filter', type = int, required = True, help = 'Min common SNP cut of value for filtering observed kinship results')
args = parser.parse_args()

if args.type != 'cons' and args.type != 'var':
	print('Please input a valid simulation type: "cons" for constant size, "var" variable size')
	exit()

if args.type == 'cons':
	if args.fam is None:
		parser.error('Family size is required for constant family size simulations!')

cut_of = str(args.filter)
iter_s = str(args.iter)
siter_s = str(args.siter)

#Functions
def random_sample():
	sim_sample = pd.DataFrame(columns = ['Period_Building', 'Period', 'Family', 'Age', 'Ind'])	
	buildings = list(data['Period_Building'])
	for b in buildings:
		build = [b]
		rand = sim[sim['Period_Building'].isin(build)].reset_index(drop = True)
		rand_sub = rand[rand['Age'] == 'Subadult']
		select_sub = rand_sub.sample(n = sim_sub[b])
		rand_adu = rand[rand['Age'] == 'Adult']
		select_adu = rand_adu.sample(n = sim_adu[b])
		rand_ado = rand[rand['Age'] == 'Adolescent']
		select_ado = rand_ado.sample(n = sim_ado[b])
		select = pd.concat([select_sub, select_adu, select_ado])
		sim_sample = pd.concat([sim_sample, select]).reset_index(drop = True)
	return sim_sample

def random_pair():
	sim_pair = pd.DataFrame(columns = ['Ind1', 'Ind2', 'Period_Building1', 'Period1', 'Family1', 'Period_Building2', 'Period2', 'Family2', 'Relatedness'])
	buildings = list(pair['Period_Building'])
	for b in buildings:
		building = [b]
		randp = sim_rand[sim_rand['Period_Building1'].isin(building)].reset_index(drop = True)

		randp_adu_adu = randp[randp['Age_Cat'] == 'Adult_Adult']
		selectp_adu_adu = randp_adu_adu.sample(n = pair_adu_adu[b])
		randp_adu_sub = randp[randp['Age_Cat'] == 'Adult_Subadult']
		selectp_adu_sub = randp_adu_sub.sample(n = pair_adu_sub[b])
		randp_sub_sub = randp[randp['Age_Cat'] == 'Subadult_Subadult']
		selectp_sub_sub = randp_sub_sub.sample(n = pair_sub_sub[b])
		randp_ado_sub = randp[randp['Age_Cat'] == 'Adolescent_Subadult']
		selectp_ado_sub = randp_ado_sub.sample(n = pair_ado_sub[b])
		randp_ado_adu = randp[randp['Age_Cat'] == 'Adolescent_Adult']
		selectp_ado_adu = randp_ado_adu.sample(n = pair_ado_adu[b])
		randp_ado_ado = randp[randp['Age_Cat'] == 'Adolescent_Adolescent']
		selectp_ado_ado = randp_ado_ado.sample(n = pair_ado_ado[b])

		selectp = pd.concat([selectp_adu_adu, selectp_adu_sub, selectp_sub_sub, selectp_ado_sub, selectp_ado_adu, selectp_ado_ado])
		sim_pair = pd.concat([sim_pair, selectp]).reset_index(drop = True)
	return sim_pair

#Read burial information for random selection
data_all = pd.read_csv('burial_DNA_info.txt', sep = '\t', converters = {'Building': str})
if not os.path.isfile('burial_DNA_info.txt'):
	raise FileNotFoundError('Burial information is not found. Please run samplecobureidkin.py first!')

data = data_all[(data_all['Period'] == 'Early') | (data_all['Period'] == 'Middle') | (data_all['Period'] == 'Late')]
data = data[data['Burials'] != 1]
data = data[data['DNA'] != 1]
data = data[data['DNA'] != 0]
data = data[data['Building'] != '']
data.dropna(subset = ['Building'], inplace = True)
data['Period_Building'] = data[['Period', 'Building']].agg('_'.join, axis=1)

sim_sub = pd.Series(data.Subadult_DNA.values, index = data.Period_Building).to_dict()
sim_adu = pd.Series(data.Adult_DNA.values, index = data.Period_Building).to_dict()
sim_ado = pd.Series(data.Adolescent_DNA.values, index = data.Period_Building).to_dict()

#Save coburial pair information for random selection
if not os.path.isfile('period_pair_info_filt' + cut_of + '.txt'):
	raise FileNotFoundError('Observed relatedness table is not found. Please run samplecobureidkin.py first or check input parameters!')
	exit()

pair = pd.read_csv('period_pair_info_filt' + cut_of + '.txt', sep = '\t', converters = {'Building': str})

pair['Period_Building'] = pair[['Period', 'Building']].agg('_'.join, axis=1)

#Get buildings with date information
dates = pd.read_csv(args.date, sep = '\t', converters = {'Building': str})

dates['Period_Building'] = dates[['Period', 'Building']].agg('_'.join, axis=1)

dates = dates[dates['Period_Building'].isin(data['Period_Building'])].reset_index(drop = True)

date_info = dates[['Period_Building', 'Date']]

building_list = date_info['Period_Building']

#Read simulated family table created
if args.type == 'cons':
	fam_s = str(args.fam)
else:
	fam_s = 'var'

if not os.path.isfile('simulation_main_distribution_fam' + fam_s +'_sim' + siter_s + '.txt'):
	raise FileNotFoundError('Simulated family table is not found. Please run samplecobureidkin.py first or check input parameters!')
	exit()

sim = pd.read_csv('simulation_main_distribution_fam' + fam_s +'_sim' + siter_s + '.txt', sep = '\t', converters = {'Building': str})

sim['Period_Building'] = sim[['Period', 'Building']].agg('_'.join, axis=1)

pair_adu_adu = pd.Series(pair.Adult_Adult.values, index = pair.Period_Building).to_dict()
pair_adu_sub = pd.Series(pair.Adult_Subadult.values, index = pair.Period_Building).to_dict()
pair_sub_sub = pd.Series(pair.Subadult_Subadult.values, index = pair.Period_Building).to_dict()
pair_ado_sub = pd.Series(pair.Adolescent_Subadult.values, index = pair.Period_Building).to_dict()
pair_ado_adu = pd.Series(pair.Adolescent_Adult.values, index = pair.Period_Building).to_dict()
pair_ado_ado = pd.Series(pair.Adolescent_Adolescent.values, index = pair.Period_Building).to_dict()

#Create empty tables to store randomization loop outputs
rand = pd.DataFrame(columns = ['Period_Building', 'Family', 'Ind', 'Age', 'Iteration'])
sim_rln = pd.DataFrame(columns = ['Ind1', 'Ind2', 'Period_Building1', 'Period1', 'Family1', 'Age1', 'Period_Building2', 'Period2', 'Family2', 'Age2', 'Relatedness'])
pair_rln = pd.DataFrame(columns = ['Ind1', 'Ind2', 'Period_Building1', 'Period1', 'Family1', 'Age1', 'Period_Building2', 'Period2', 'Family2', 'Age2', 'Relatedness'])

rand_freq = pd.DataFrame(columns = ['Period', 'Related', 'Unrelated', 'Iteration'])
rand_freq_sub = pd.DataFrame(columns = ['Period', 'Related', 'Unrelated', 'Iteration'])
rand_freq_sub_ado = pd.DataFrame(columns = ['Period', 'Related', 'Unrelated', 'Iteration'])

randb_freq = pd.DataFrame(columns = ['Period_Building', 'Related', 'Unrelated', 'Iteration'])
randb_freq_sub = pd.DataFrame(columns = ['Period_Building', 'Related', 'Unrelated', 'Iteration'])
randb_freq_sub_ado = pd.DataFrame(columns = ['Period_Building', 'Related', 'Unrelated', 'Iteration'])

#Start randomization
for i in range(1, (args.iter + 1)):
	sample_iter = random_sample()
	i_s = str(i)
	sample_iter['Iteration'] = i_s
	rand = pd.concat([rand, sample_iter])

	sample_rel = sample_iter.copy()
	inds_rand = sample_rel['Ind'].to_list()
	comb_rand = pd.DataFrame(list(product(inds_rand, inds_rand)))
	comb_rand[2] = ['-'.join(x) for x in np.sort(comb_rand[[0, 1]], axis = 1)]
	comb_rand = comb_rand.drop_duplicates(subset = [2], keep = 'first')
	comb_rand = comb_rand.drop([2], axis = 1);
	ind1_rand = pd.merge(comb_rand, sample_rel, left_on = 0, right_on = 'Ind', how = 'inner')
	ind1_rand = ind1_rand.rename(columns={0: 'Ind1'})
	sim_rand = pd.merge(ind1_rand, sample_rel, left_on = 1, right_on = 'Ind', how = 'inner')
	sim_rand = sim_rand.rename(columns={1: 'Ind2'})
	sim_rand = sim_rand.drop(['Ind_x', 'Ind_y'], axis = 1)
	sim_rand.columns = sim_rand.columns.str.replace('_x', '1', regex = True)
	sim_rand.columns = sim_rand.columns.str.replace('_y', '2', regex = True)
	sim_rand = sim_rand[sim_rand['Period1'] == sim_rand['Period2']]
	sim_rand['Relatedness'] = ''
	sim_rand.loc[sim_rand['Family1'] == sim_rand['Family2'], 'Relatedness'] = 'Related'
	sim_rand.loc[sim_rand['Family1'] != sim_rand['Family2'], 'Relatedness'] = 'Unrelated'
	sim_rand = sim_rand.drop(['Iteration1'], axis = 1)
	sim_rand = sim_rand.rename(columns={'Iteration2': 'Iteration'})
	sim_rand = sim_rand[sim_rand['Ind1'] != sim_rand['Ind2']]
	sim_rln = pd.concat([sim_rln, sim_rand])

	sim_rand['Age_Cat'] = ['_'.join(x) for x in np.sort(sim_rand[['Age1', 'Age2']], axis=1)]
	sim_rand = sim_rand[(sim_rand['Period_Building1'] == sim_rand['Period_Building2'])]
	pair_iter = random_pair()
	pair_rln = pd.concat([pair_rln, pair_iter])

	pair_iter.loc[(pair_iter.Age1 == 'Subadult'), 'Age1_Cat_Ado_Adult'] = 'Subadult'
	pair_iter.loc[(pair_iter.Age1 == 'Adult'), 'Age1_Cat_Ado_Adult'] = 'Adult'
	pair_iter.loc[(pair_iter.Age1 == 'Adolescent'), 'Age1_Cat_Ado_Adult'] = 'Adult'
	pair_iter.loc[(pair_iter.Age2 == 'Subadult'), 'Age2_Cat_Ado_Adult'] = 'Subadult'
	pair_iter.loc[(pair_iter.Age2 == 'Adult', 'Age2_Cat_Ado_Adult')] = 'Adult'
	pair_iter.loc[(pair_iter.Age2 == 'Adolescent', 'Age2_Cat_Ado_Adult')] = 'Adult'

	pair_iter.loc[(pair_iter.Age1 == 'Subadult', 'Age1_Cat_Ado_Subadult')] = 'Subadult'
	pair_iter.loc[(pair_iter.Age1 == 'Adult', 'Age1_Cat_Ado_Subadult')] = 'Adult'
	pair_iter.loc[(pair_iter.Age1 == 'Adolescent', 'Age1_Cat_Ado_Subadult')] = 'Subadult'
	pair_iter.loc[(pair_iter.Age2 == 'Subadult', 'Age2_Cat_Ado_Subadult')] = 'Subadult'
	pair_iter.loc[(pair_iter.Age2 == 'Adult', 'Age2_Cat_Ado_Subadult')] = 'Adult'
	pair_iter.loc[(pair_iter.Age2 == 'Adolescent', 'Age2_Cat_Ado_Subadult')] = 'Subadult'

	freq_table = pair_iter[['Period1', 'Relatedness']]
	freq_table = freq_table.rename(columns={'Period1': 'Period'})
	freq_iter = freq_table.pivot_table(index = 'Period', columns = 'Relatedness', aggfunc = lambda x: len(x))
	freq_iter = freq_iter.reindex(['Early', 'Middle', 'Late'])
	freq_iter = freq_iter.fillna(0)
	freq_iter = freq_iter.astype(int)
	freq_iter = freq_iter.reset_index().rename_axis(None, axis = 1)
	freq_iter['Iteration'] = i_s
	rand_freq = pd.concat([rand_freq, freq_iter])

	freq_table_sub = pair_iter[(pair_iter['Age1_Cat_Ado_Adult'] == pair_iter['Age2_Cat_Ado_Adult'])]
	freq_table_sub = freq_table_sub[freq_table_sub['Age1_Cat_Ado_Adult'] == 'Subadult']
	freq_table_sub = freq_table_sub[['Period1', 'Relatedness']]
	freq_table_sub = freq_table_sub.rename(columns={'Period1': 'Period'})
	freq_iter_sub = freq_table_sub.pivot_table(index = 'Period', columns = 'Relatedness', aggfunc = lambda x: len(x))
	freq_iter_sub = freq_iter_sub.reindex(['Early', 'Middle', 'Late'])
	freq_iter_sub = freq_iter_sub.fillna(0)
	freq_iter_sub = freq_iter_sub.astype(int)
	freq_iter_sub = freq_iter_sub.reset_index().rename_axis(None, axis = 1)
	freq_iter_sub['Iteration'] = i_s
	rand_freq_sub = pd.concat([rand_freq_sub, freq_iter_sub])

	freq_table_sub_ado = pair_iter[(pair_iter['Age1_Cat_Ado_Subadult'] == pair_iter['Age2_Cat_Ado_Subadult'])]
	freq_table_sub_ado = freq_table_sub_ado[freq_table_sub_ado['Age1_Cat_Ado_Subadult'] == 'Subadult'] 
	freq_table_sub_ado = freq_table_sub_ado[['Period1', 'Relatedness']]
	freq_table_sub_ado = freq_table_sub_ado.rename(columns={'Period1': 'Period'})
	freq_iter_sub_ado = freq_table_sub_ado.pivot_table(index = 'Period', columns = 'Relatedness', aggfunc = lambda x: len(x))
	freq_iter_sub_ado = freq_iter_sub_ado.reindex(['Early', 'Middle', 'Late'])
	freq_iter_sub_ado = freq_iter_sub_ado.fillna(0)
	freq_iter_sub_ado = freq_iter_sub_ado.astype(int)
	freq_iter_sub_ado = freq_iter_sub_ado.reset_index().rename_axis(None, axis = 1)
	freq_iter_sub_ado['Iteration'] = i_s
	rand_freq_sub_ado = pd.concat([rand_freq_sub_ado, freq_iter_sub_ado])

	freqb_table = pair_iter[['Period_Building1', 'Relatedness']]
	freqb_table = freqb_table.rename(columns={'Period_Building1': 'Period_Building'})
	freqb_iter = freqb_table.pivot_table(index = 'Period_Building', columns = 'Relatedness', aggfunc = lambda x: len(x))
	freqb_iter = freqb_iter.reset_index().rename_axis(None, axis = 1)
	freqb_list = freqb_iter['Period_Building']
	diff = list(set(building_list) ^ set(freqb_list))
	if diff != 0:
		for i in diff:
			values = [i, 0, 0,]
			df = pd.DataFrame(columns=['Period_Building', 'Related', 'Unrelated'], data=[values])
			freqb_iter = pd.concat([freqb_iter, df], axis = 0).reset_index(drop = True)
	freqb_iter = freqb_iter.fillna(0).set_index('Period_Building').reindex(building_list).reset_index()
	freqb_iter['Iteration'] = i_s
	randb_freq = pd.concat([randb_freq, freqb_iter])
	
	freqb_table_sub = pair_iter[(pair_iter['Age1_Cat_Ado_Adult'] == pair_iter['Age2_Cat_Ado_Adult'])]
	freqb_table_sub = freqb_table_sub[freqb_table_sub['Age1_Cat_Ado_Adult'] == 'Subadult']
	freqb_table_sub = freqb_table_sub[['Period_Building1', 'Relatedness']]
	freqb_table_sub = freqb_table_sub.rename(columns={'Period_Building1': 'Period_Building'})
	freqb_iter_sub = freqb_table_sub.pivot_table(index = 'Period_Building', columns = 'Relatedness', aggfunc = lambda x: len(x))
	freqb_iter_sub = freqb_iter_sub.reset_index().rename_axis(None, axis = 1)
	freqb_list = freqb_iter_sub['Period_Building']
	diff = list(set(building_list) ^ set(freqb_list))
	if diff != 0:
		for i in diff:
			values = [i, 0, 0,]
			df = pd.DataFrame(columns=['Period_Building', 'Related', 'Unrelated'], data=[values])
			freqb_iter_sub = pd.concat([freqb_iter_sub, df], axis = 0).reset_index(drop = True)
	freqb_iter_sub = freqb_iter_sub.fillna(0).set_index('Period_Building').reindex(building_list).reset_index()
	freqb_iter_sub['Iteration'] = i_s
	randb_freq_sub = pd.concat([randb_freq_sub, freqb_iter_sub])

	freqb_table_sub_ado = pair_iter[(pair_iter['Age1_Cat_Ado_Subadult'] == pair_iter['Age2_Cat_Ado_Subadult'])]
	freqb_table_sub_ado = freqb_table_sub_ado[freqb_table_sub_ado['Age1_Cat_Ado_Subadult'] == 'Subadult'] 
	freqb_table_sub_ado = freqb_table_sub_ado[['Period_Building1', 'Relatedness']]
	freqb_table_sub_ado = freqb_table_sub_ado.rename(columns={'Period_Building1': 'Period_Building'})
	freqb_iter_sub_ado = freqb_table_sub_ado.pivot_table(index = 'Period_Building', columns = 'Relatedness', aggfunc = lambda x: len(x))
	freqb_iter_sub_ado = freqb_iter_sub_ado.reset_index().rename_axis(None, axis = 1)
	freqb_list = freqb_iter_sub_ado['Period_Building']
	diff = list(set(building_list) ^ set(freqb_list))
	if diff != 0:
		for i in diff:
			values = [i, 0, 0,]
			df = pd.DataFrame(columns=['Period_Building', 'Related', 'Unrelated'], data=[values])
			freqb_iter_sub_ado = pd.concat([freqb_iter_sub_ado, df], axis = 0).reset_index(drop = True)
	freqb_iter_sub_ado = freqb_iter_sub_ado.fillna(0).set_index('Period_Building').reindex(building_list).reset_index()
	freqb_iter_sub_ado['Iteration'] = i_s
	randb_freq_sub_ado = pd.concat([randb_freq_sub_ado, freqb_iter_sub_ado])

#Organize ands save output tables
#List of individuals selected after randomization
rand['Iteration'] = rand['Iteration'].astype(str)
rand[['Period', 'Building']] = rand.Period_Building.str.split('_', expand = True)
rand = rand.drop(['Period_Building'], axis = 1)
rand['Simulation'] = siter_s
rand = rand.iloc[:, [6,3,4,5,1,0,2]]
rand.to_csv('simulation_main_distribution_fam' + fam_s + '_sim' + siter_s + '_iter' + iter_s + '_filt' + cut_of + '.txt', sep = '\t', index = False)

#Relatedness table of individuals selected after randomization
sim_rln['Iteration'] = sim_rln['Iteration'].astype(str)
sim_rln = sim_rln.drop(['Period1', 'Period2'], axis = 1)
sim_rln[['Period1', 'Building1']] = sim_rln.Period_Building1.str.split('_', expand = True)
sim_rln[['Period2', 'Building2']] = sim_rln.Period_Building2.str.split('_', expand = True)
sim_rln = sim_rln.drop(['Period_Building1', 'Period_Building2'], axis = 1)
sim_rln['Simulation'] = siter_s
sim_rln = sim_rln.iloc[:, [12,9,10,7,0,2,3,11,8,1,4,5,6]]
sim_rln.to_csv('simulation_sample_relatedness_fam' + fam_s + '_sim' + siter_s + '_iter' + iter_s + '_filt' + cut_of + '.txt', sep = '\t', index = False)

#List of pairs selected after randomization
pair_rln['Iteration'] = pair_rln['Iteration'].astype(str)
pair_rln = pair_rln.drop(['Period1', 'Period2'], axis = 1)
pair_rln[['Period1', 'Building1']] = pair_rln.Period_Building1.str.split('_', expand = True)
pair_rln[['Period2', 'Building2']] = pair_rln.Period_Building2.str.split('_', expand = True)
pair_rln = pair_rln.drop(['Period_Building1', 'Period_Building2', 'Age_Cat'], axis = 1)
pair_rln['Simulation'] = siter_s
pair_rln = pair_rln.iloc[:, [12,9,10,7,0,2,3,11,8,1,4,5,6]]
pair_rln.to_csv('simulation_pair_relatedness_fam' + fam_s + '_sim' + siter_s + '_iter' + iter_s + '_filt' + cut_of + '.txt', sep = '\t', index = False)

#Total number of related and unrelated pairs per period for all samples
rand_freq['Iteration'] = rand_freq['Iteration'].astype('str')
rand_freq['Total'] = rand_freq['Related'] + rand_freq['Unrelated']
rand_freq['Simulation'] = siter_s
rand_freq = rand_freq.iloc[:, [5,3,0,1,2,4]].replace('',0).fillna(0)
rand_freq['Related'] = rand_freq['Related'].astype(int)
rand_freq['Unrelated'] = rand_freq['Unrelated'].astype(int)
rand_freq['Total'] = rand_freq['Total'].astype(int)
rand_freq.to_csv('simulation_random_relatedness_freq_fam' + fam_s + '_sim' + siter_s + '_iter' + iter_s + '_filt' + cut_of + '.txt', sep = '\t', index = False)

#Total number of related and unrelated pairs per period for subadults exc adolescents
rand_freq_sub['Iteration'] = rand_freq_sub['Iteration'].astype('str')
rand_freq_sub['Total'] = rand_freq_sub['Related'] + rand_freq_sub['Unrelated']
rand_freq_sub['Simulation'] = siter_s
rand_freq_sub = rand_freq_sub.iloc[:, [5,3,0,1,2,4]].replace('',0).fillna(0)
rand_freq_sub['Related'] = rand_freq_sub['Related'].astype(int)
rand_freq_sub['Unrelated'] = rand_freq_sub['Unrelated'].astype(int)
rand_freq_sub['Total'] = rand_freq_sub['Total'].astype(int)
rand_freq_sub.to_csv('simulation_random_relatedness_freq_subadults_wo_ados_fam' + fam_s + '_sim' + siter_s + '_iter' + iter_s + '_filt' + cut_of + '.txt', sep = '\t', index = False)

#Total number of related and unrelated pairs per period for subadults inc adolescents
rand_freq_sub_ado['Iteration'] = rand_freq_sub_ado['Iteration'].astype('str')
rand_freq_sub_ado['Total'] = rand_freq_sub_ado['Related'] + rand_freq_sub_ado['Unrelated']
rand_freq_sub_ado['Simulation'] = siter_s
rand_freq_sub_ado = rand_freq_sub_ado.iloc[:, [5,3,0,1,2,4]].replace('',0).fillna(0)
rand_freq_sub_ado['Related'] = rand_freq_sub_ado['Related'].astype(int)
rand_freq_sub_ado['Unrelated'] = rand_freq_sub_ado['Unrelated'].astype(int)
rand_freq_sub_ado['Total'] = rand_freq_sub_ado['Total'].astype(int)
rand_freq_sub_ado.to_csv('simulation_random_relatedness_freq_subadults_wt_ados_fam' + fam_s + '_sim' + siter_s + '_iter' + iter_s + '_filt' + cut_of + '.txt', sep = '\t', index = False)

randb_freq['Iteration'] = randb_freq['Iteration'].astype('str')
randb_freq['Total'] = randb_freq['Related'] + randb_freq['Unrelated']
randb_freq['Simulation'] = siter_s
randb_freq = randb_freq.merge(date_info, on = 'Period_Building', how = 'left')
randb_freq[['Period', 'Building']] = randb_freq.Period_Building.str.split('_', expand = True)
randb_freq = randb_freq.drop(['Period_Building'], axis = 1)
randb_freq = randb_freq.iloc[:, [4,2,6,7,5,0,1,3]].replace('',0).fillna(0)
randb_freq['Related'] = randb_freq['Related'].astype(int)
randb_freq['Unrelated'] = randb_freq['Unrelated'].astype(int)
randb_freq['Total'] = randb_freq['Total'].astype(int)
randb_freq.to_csv('simulation_random_relatedness_per_building_freq_fam' + fam_s + '_sim' + siter_s + '_iter' + iter_s + '_filt' + cut_of + '.txt', sep = '\t', index = False)

randb_freq_sub['Iteration'] = randb_freq_sub['Iteration'].astype('str')
randb_freq_sub['Total'] = randb_freq_sub['Related'] + randb_freq_sub['Unrelated']
randb_freq_sub['Simulation'] = siter_s
randb_freq_sub = randb_freq_sub.merge(date_info, on = 'Period_Building', how = 'left')
randb_freq_sub[['Period', 'Building']] = randb_freq_sub.Period_Building.str.split('_', expand = True)
randb_freq_sub = randb_freq_sub.drop(['Period_Building'], axis = 1)
randb_freq_sub = randb_freq_sub.iloc[:, [4,2,6,7,5,0,1,3]].replace('',0).fillna(0)
randb_freq_sub['Related'] = randb_freq_sub['Related'].astype(int)
randb_freq_sub['Unrelated'] = randb_freq_sub['Unrelated'].astype(int)
randb_freq_sub['Total'] = randb_freq_sub['Total'].astype(int)
randb_freq_sub.to_csv('simulation_random_relatedness_per_building_freq_subadults_wo_ados_fam' + fam_s + '_sim' + siter_s + '_iter' + iter_s + '_filt' + cut_of + '.txt', sep = '\t', index = False)

randb_freq_sub_ado['Iteration'] = randb_freq_sub_ado['Iteration'].astype('str')
randb_freq_sub_ado['Total'] = randb_freq_sub_ado['Related'] + randb_freq_sub_ado['Unrelated']
randb_freq_sub_ado['Simulation'] = siter_s
randb_freq_sub_ado = randb_freq_sub_ado.merge(date_info, on = 'Period_Building', how = 'left')
randb_freq_sub_ado[['Period', 'Building']] = randb_freq_sub_ado.Period_Building.str.split('_', expand = True)
randb_freq_sub_ado = randb_freq_sub_ado.drop(['Period_Building'], axis = 1)
randb_freq_sub_ado = randb_freq_sub_ado.iloc[:, [4,2,6,7,5,0,1,3]].replace('',0).fillna(0)
randb_freq_sub_ado['Related'] = randb_freq_sub_ado['Related'].astype(int)
randb_freq_sub_ado['Unrelated'] = randb_freq_sub_ado['Unrelated'].astype(int)
randb_freq_sub_ado['Total'] = randb_freq_sub_ado['Total'].astype(int)
randb_freq_sub_ado.to_csv('simulation_random_relatedness_per_building_freq_subadults_wt_ados_fam' + fam_s + '_sim' + siter_s + '_iter' + iter_s + '_filt' + cut_of + '.txt', sep = '\t', index = False)

print('End of random selection!!')
