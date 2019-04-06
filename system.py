import numpy as np
from Bio.Seq import Seq
from scipy import stats
from sklearn import datasets, linear_model
from scipy.special import binom
from time import time
import math

def Import_Data(File):
	
	global codons, amino_acids, freqs, amino_acid_codons, N, organism, file_name, num_data, rev_comp, trans
	
	file_name = 'Data/'+File+'.data'
	organism = file_name.split('/')[-1].split('_ppm')[0]
	
	freqs = {}
	codons = set([])
	amino_acids = set([])
	data_type = ''
	num_data = 0
	
	in_file = open(file_name, 'r')
	
	for line in in_file.readlines():
		
		if not line[0] == '#':
			term = line.split(':')[0]
			
			if len(term) == 3:
				codons.add(term)
			
			if len(term) < 3:
				amino_acids.add(term)
		
			key = data_type + term
			
			if key in freqs:
				freqs[key] += float(line.split(':')[1])
			else:
				freqs[key] = float(line.split(':')[1])
		
		else: 
			if len(line.split()) == 7 and line.split()[4] == 'set':
				data_type = line.split()[5][:-1]+':'
				num_data += 1/2
			else:
				data_type = ''
	in_file.close()
	
	codons = list(codons)
	codons.sort()
	
	rev_comp = {}
	trans = {}
	
	for codon in codons:
		rev_comp[codon] = str(Seq(codon).reverse_complement())
		trans[codon] = str(Seq(codon).translate())
		
	amino_acids = list(amino_acids)
	amino_acid_codons = {amino_acid:[] for amino_acid in amino_acids}
	
	for codon in codons:
		amino_acid = str(Seq(codon).translate())
	
		if amino_acid in amino_acids:
			amino_acid_codons[amino_acid].append(codon)

	amino_acids = list(set(amino_acids))
	amino_acids.sort()

	
	N = len(codons)
	
	
def Set_Data(set_num):
	
	global codon_freqs, amino_freqs, freqs, data_points, best_logL, data_set
	
	data_set = set_num
	
	if data_set == 'all':
		codon_freqs = np.zeros(N)
		amino_freqs = np.zeros(len(amino_acids))
		
		for i in range(int(num_data)):
			codon_freqs += np.array([freqs[str(i) + ':' + codon] for codon in codons])
			amino_freqs += np.array([freqs[str(i) + ':' + amino] for amino in amino_acids])/num_data
			
	else:
		codon_freqs = np.array([freqs[str(data_set) + ':' + codon] for codon in codons])
		amino_freqs = np.array([freqs[str(data_set) + ':' + amino] for amino in amino_acids])
	
	data_points = sum(codon_freqs)
	codon_freqs = codon_freqs/data_points

	
def Generate_Seq():
	
	if n_params == 3:
		return [str(round(init_ranges[i]*term,2)) for i,term in enumerate(np.random.random(n_params))] + ['0']*2
	else:
		return [str(round(init_ranges[i]*term,2)) for i,term in enumerate(np.random.random(n_params))]

dx = 1e2
depth = 2
curr_depth = -1

def Set_Params(model_params):
	global kappa1_i, kappa2_i, sbeta_i, T0beta_i, r_i, mutables, init_ranges, max_vals, min_vals, n_params
	n_params = model_params
	
	if model_params == 3:
		r_i = [4, 5]
		
	elif model_params == 5:

		r_i = [4, 5]
		
	elif model_params == 7:

		r_i = [4, 7]
		
	elif model_params == 12:

		r_i = [4, 12]
		
	elif model_params == 16:

		r_i = [4, 16]
		
	elif model_params == 19:

		r_i = [3, model_params]
	kappa1_i = 0
	kappa2_i = 1
	sbeta_i = 2
	
	if not model_params == 19:
		T0beta_i = 3
		
	mutables = [i for i in range(model_params)]
	init_ranges = [100 for i in range(r_i[-1])]
	max_vals = [10**2 for i in range(r_i[-1])]
	min_vals = [0.1, 0.1] + [0.0 for i in range(2, r_i[-1])]

def is_Transition(codon_1, codon_2):
	
	trans = False
	del_nuc = ''
	
	for k in range(len(codon_1)):
		if not codon_1[k] == codon_2[k]:
			del_nuc = codon_1[k]
			if (codon_1[k] in ['T', 'C'] and codon_2[k] in ['T', 'C']) or\
				(codon_1[k] in ['G', 'A'] and codon_2[k] in ['G', 'A']):
				trans = True
	return trans, del_nuc
		
def Get_Ceff(i,rs):
	
	if len(rs) == 1:
		
		r =	{'A/A': rs[0],\
			 'A/C': rs[0],\
			 'A/G': rs[0],\
			 'A/T': 1.,\

			 'C/A': rs[0],\
			 'C/C': 0.,\
			 'C/G': 1.,\
			 'C/T': 0.,\

			 'G/A': 0.,\
			 'G/C': 1.,\
			 'G/G': 0.,\
			 'G/T': rs[0],\

			 'T/A': 1.,\
			 'T/C': rs[0],\
			 'T/G': rs[0],\
			 'T/T': rs[0]}
	
	elif len(rs) == 3:
		
		r =	{'A/A': rs[1],\
			 'A/C': rs[0],\
			 'A/G': rs[2],\
			 'A/T': 1.,\

			 'C/A': rs[0],\
			 'C/C': 0.,\
			 'C/G': 1.,\
			 'C/T': 0.,\

			 'G/A': 0.,\
			 'G/C': 1.,\
			 'G/G': 0.,\
			 'G/T': rs[0],\

			 'T/A': 1.,\
			 'T/C': rs[2],\
			 'T/G': rs[0],\
			 'T/T': rs[1]}
			 
	elif len(rs) == 8:
		
		r =	{'A/A': rs[0],\
			 'A/C': rs[1],\
			 'A/G': rs[2],\
			 'A/T': 1.,\

			 'C/A': rs[3],\
			 'C/C': 1e-9,\
			 'C/G': 1.,\
			 'C/T': 1e-9,\

			 'G/A': 1e-9,\
			 'G/C': 1.,\
			 'G/G': 1e-9,\
			 'G/T': rs[4],\

			 'T/A': 1.,\
			 'T/C': rs[5],\
			 'T/G': rs[6],\
			 'T/T': rs[7]}
	
	elif len(rs) == 12:
		
		r =	{'A/A': rs[0],\
			 'A/C': rs[1],\
			 'A/G': rs[2],\
			 'A/T': 1.,\

			 'C/A': rs[3],\
			 'C/C': rs[4],\
			 'C/G': 1.,\
			 'C/T': rs[5],\

			 'G/A': rs[6],\
			 'G/C': 1.,\
			 'G/G': rs[7],\
			 'G/T': rs[8],\

			 'T/A': 1.,\
			 'T/C': rs[9],\
			 'T/G': rs[10],\
			 'T/T': rs[11]}
	
	elif len(rs) == 16:
		
		r =	{'A/A': rs[0],\
			 'A/C': rs[1],\
			 'A/G': rs[2],\
			 'A/T': rs[3],\

			 'C/A': rs[4],\
			 'C/C': rs[5],\
			 'C/G': rs[6],\
			 'C/T': rs[7],\

			 'G/A': rs[8],\
			 'G/C': rs[9],\
			 'G/G': rs[10],\
			 'G/T': rs[11],\

			 'T/A': rs[12],\
			 'T/C': rs[13],\
			 'T/G': rs[14],\
			 'T/T': rs[15]}
	
	codon = codons[i]
	cog_nuc = rev_comp[codon][0]
	anticodon23 = rev_comp[codon][1:]
	
	C_eff = 0
	
	for nuc in ['A', 'C', 'G', 'T']:
	
		anticodon = nuc + anticodon23
		
		C_eff += r[nuc + '/' + codon[2]]*freqs['tr'+anticodon]
	
	return C_eff
	
def Get_sCeff(i,rs,A,s):
	
	if len(rs) == 1:
		
		r =	{'A/A': rs[0],\
			 'A/C': rs[0],\
			 'A/G': rs[0],\
			 'A/T': 1.,\

			 'C/A': rs[0],\
			 'C/C': 0.,\
			 'C/G': 1.,\
			 'C/T': 0.,\

			 'G/A': 0.,\
			 'G/C': 1.,\
			 'G/G': 0.,\
			 'G/T': rs[0],\

			 'T/A': 1.,\
			 'T/C': rs[0],\
			 'T/G': rs[0],\
			 'T/T': rs[0]}
	
	elif len(rs) == 3:
		
		r =	{'A/A': rs[1],\
			 'A/C': rs[0],\
			 'A/G': rs[2],\
			 'A/T': 1.,\

			 'C/A': rs[0],\
			 'C/C': 0.,\
			 'C/G': 1.,\
			 'C/T': 0.,\

			 'G/A': 0.,\
			 'G/C': 1.,\
			 'G/G': 0.,\
			 'G/T': rs[0],\

			 'T/A': 1.,\
			 'T/C': rs[2],\
			 'T/G': rs[0],\
			 'T/T': rs[1]}
			 
	elif len(rs) == 8:
		
		r =	{'A/A': rs[0],\
			 'A/C': rs[1],\
			 'A/G': rs[2],\
			 'A/T': 1.,\

			 'C/A': rs[3],\
			 'C/C': 1e-9,\
			 'C/G': 1.,\
			 'C/T': 1e-9,\

			 'G/A': 1e-9,\
			 'G/C': 1.,\
			 'G/G': 1e-9,\
			 'G/T': rs[4],\

			 'T/A': 1.,\
			 'T/C': rs[5],\
			 'T/G': rs[6],\
			 'T/T': rs[7]}
	
	elif len(rs) == 12:
		
		r =	{'A/A': rs[0],\
			 'A/C': rs[1],\
			 'A/G': rs[2],\
			 'A/T': 1.,\

			 'C/A': rs[3],\
			 'C/C': rs[4],\
			 'C/G': 1.,\
			 'C/T': rs[5],\

			 'G/A': rs[6],\
			 'G/C': 1.,\
			 'G/G': rs[7],\
			 'G/T': rs[8],\

			 'T/A': 1.,\
			 'T/C': rs[9],\
			 'T/G': rs[10],\
			 'T/T': rs[11]}
	
	elif len(rs) == 16:
		
		r =	{'A/A': rs[0],\
			 'A/C': rs[1],\
			 'A/G': rs[2],\
			 'A/T': rs[3],\

			 'C/A': rs[4],\
			 'C/C': rs[5],\
			 'C/G': rs[6],\
			 'C/T': rs[7],\

			 'G/A': rs[8],\
			 'G/C': rs[9],\
			 'G/G': rs[10],\
			 'G/T': rs[11],\

			 'T/A': rs[12],\
			 'T/C': rs[13],\
			 'T/G': rs[14],\
			 'T/T': rs[15]}
	
	codon = codons[i]
	
	anticodon23 = str(rev_comp[codon])[1:]
	
	sC_eff = 0
	
	for nuc in ['A', 'C', 'G', 'T']:
	
		anticodon = nuc + anticodon23
		
		cog_A = str(trans[rev_comp[anticodon]])
			
		sC_eff += s*r[nuc + '/' + codon[2]]*freqs['tr'+anticodon]*int(not cog_A == A)
	
	return sC_eff

def Get_s(i,A,s):
	if trans[codons[i]] == A:
		return 0
	else:
		return s

def Get_Freqs(seq):
	global mut_mat, codon_freqs
	
	beta = 1e-7
	kappa1 = float(seq[kappa1_i])
	kappa2 = float(seq[kappa2_i])
	
	s = float(seq[sbeta_i])
	
	if not n_params == 19:
		T0 = float(seq[T0beta_i])*beta
	else:
		T0 = beta
		
	rs = [float(r) for r in seq[r_i[0]:r_i[1]]]
	
	pred_codon_freqs = np.zeros(N)
	
	mut_mat = np.array([[0.0 for i in range(N)] for j in range(N)])
	
	for i,codon_1 in enumerate(codons):
		for j,codon_2 in enumerate(codons):
			# codon_1 -> codon_2 / i -> j
			mut_mat[j][i] = beta*freqs['pi'+codon_2]*int(sum([int(codon_1[k] == codon_2[k]) for k in range(len(codon_1))]) == 2)
			is_tran,del_nuc = is_Transition(codon_1, codon_2)
	
			if is_tran:
				if del_nuc in ['T', 'C']:
					mut_mat[j][i] = kappa1*mut_mat[j][i]
			
				elif del_nuc in ['A', 'G']:
					mut_mat[j][i] = kappa2*mut_mat[j][i]
			
	
	for amino_acid in amino_acids:
		
		fit_mat = np.array([[0.0 for i in range(N)] for j in range(N)])
		
		for i in range(N):
			
			Ceff = Get_Ceff(i,rs)
			
			if abs(Ceff) > 0.:
				sCeff = Get_sCeff(i, rs, amino_acid, s)
				seff = sCeff/Ceff
				
				fit_mat[i][i] = 1. - seff
				fit_mat[i][i] = fit_mat[i][i]*(1. - T0/Ceff)
				
			else:
				if T0 > 0.:
					fit_mat[i][i] = 0.
				else:
					s = Get_s(i, amino_acid, s)
					fit_mat[i][i] = 1. - s
			
			
		U = np.dot(np.identity(N) + np.matrix(mut_mat), np.matrix(fit_mat))
		t1 = time()
		eigen = np.linalg.eig(U)
		
		eig_vals = list(eigen[0])
		eig_vecs = eigen[1].T
		max_eig = max(eig_vals)
		eig_index = eig_vals.index(max_eig)
		max_vec = np.array(eig_vecs[eig_index])[0]
		ps = max_vec/sum(max_vec)
		
		pred_codon_freqs = pred_codon_freqs + amino_freqs[amino_acids.index(amino_acid)]*ps
		
	return list(pred_codon_freqs)
	
def Get_Fitness(seq):
	pred_codon_freqs = Get_Freqs(seq)
	
	return -sum(abs(codon_freqs - pred_codon_freqs))/2
	
def Single_Moves(moves):
	
	nns = set([])

	for i in mutables:
		nns.add('_'.join(moves[:i] + [str(int(moves[i]) + 1)] + moves[i+1:]))
		nns.add('_'.join(moves[:i] + [str(int(moves[i]) - 1)] + moves[i+1:]))
		
	return nns

def Get_All_Neighbors(seq):
	global curr_depth, ns
	
	if not curr_depth == depth:
		curr_depth = depth
		ns = Single_Moves(['0' for i in range(len(seq))])
	
		for d in range(depth - 1):
	
			curr_ns = set(ns)
		
			for n in curr_ns:
		
				new_ns = Single_Moves(n.split('_'))
				ns = ns.union(new_ns)
		curr_ns = set(ns)
	
		for n in curr_ns:
			move_sum = sum([abs(int(move)) for move in n.split('_')])
			if move_sum < depth:
				ns.remove(n)
	
	ns = list(ns)
	ns.sort()
	print(len(seq))
	print(len(min_vals))
	print(len(max_vals))
	neighbors = [[[str(min([max([float(seq[i]) + dx*int(move), min_vals[i]]), max_vals[i]])) for i,move in enumerate(n.split('_'))],1] for n in ns]
	
	return neighbors
	

