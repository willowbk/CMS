import numpy as np
from Bio.SubsMat.MatrixInfo import blosum62 as blos
from Bio.Align.Applications import MuscleCommandline
from Bio import pairwise2
from io import StringIO
from Bio import AlignIO
from Bio.Seq import Seq

print('Importing gene sequences ... ')
import codon_genes_e_coli as gen_abun
print('Done.')
gen_abun.genes = ['PAWHEAE', 'HEAGAWGHEE']
n = len(gen_abun.genes)
print('Total genes: '+str(n))
print()

d = 8
e = 0

def s(Ai, Aj):
	#if Ai == Aj:
	#	return 2
	#else:
	#	return -1
	
	if (Ai, Aj) in blos:
		return blos[(Ai, Aj)]
	else:
		return blos[(Aj, Ai)]
		
		
		
def Pairwise_Align(seq1, seq2, align_type):
	F = np.zeros((len(seq1) + 1, len(seq2) + 1))
	dirs = np.zeros((len(seq1) + 1, len(seq2) + 1))
	max_F_i_j = [0,0,0]
	
	for i in range(int(align_type == 'local'), len(seq1) + 1):
		for j in range(int(align_type == 'local'), len(seq2) + 1):
			
			if i == 0 or j == 0:
				
				F[i][j] = -d*(i + j)
				dirs[i][j] = 2*int(i>0) + int(j>0)
				
			else:
			
				terms = []
				
				if align_type == 'local':
					terms.append([0, -3])
					
				terms.append([F[i-1][j-1] + s(seq1[i-1], seq2[j-1]), -1])
				terms.append([F[i][j-1] - d, 0])
				terms.append([F[i-1][j] - d, -2])
			
				sorted_terms = sorted(terms)
			
				F[i][j] = sorted_terms[-1][0]
				dirs[i][j] = terms.index(sorted_terms[-1]) - int(align_type == 'local')
				print(sorted_terms)
				if F[i][j] > max_F_i_j[0]:
					max_F_i_j[0] = F[i][j]
					max_F_i_j[1] = i
					max_F_i_j[2] = j
	
	print(F)
	print(dirs)
	ali1 = ''
	ali2 = ''
	score = F[i][j]

	if align_type == 'local':
		score = max_F_i_j[0]
		i = max_F_i_j[1]
		j = max_F_i_j[2]	
		print()
	
	while i + j > 0:
	
		if dirs[i][j] == 0:
			ali1 = seq1[i-1] + ali1 
			ali2 = seq2[j-1] + ali2
		
			i -= 1
			j -= 1
		
		elif dirs[i][j] == 1:
			ali1 = '-' + ali1
			ali2 = seq2[j-1] + ali2
		
			j -= 1
	
		elif dirs[i][j] == 2:
			ali1 = seq1[i-1] + ali1
			ali2 = '-' + ali2 
		
			i -= 1
		
		elif dirs[i][j] == -1:
			i = 0
			j = 0
		
	return ali1,ali2,score

for i in range(n):
	for j in range(i + 1, n):
		seq1 = gen_abun.genes[i]
		seq2 = gen_abun.genes[j]
		
		disp = 40
		
		print('\tgene '+str(i)+': '+seq1[:disp] + '...'*int(len(seq1) > disp))
		print('\tgene '+str(j)+': '+seq2[:disp] + '...'*int(len(seq2) > disp))
		print()
		ali1,ali2,score = Pairwise_Align(seq1, seq2, 'global')
		print('\tAlignment:')
		print('\t'+ali1[:disp] + '...'*int(len(ali1) > disp)) 
		print('\t'+ali2[:disp] + '...'*int(len(ali2) > disp)) 
		print('\tscore = '+str(score))
		print()
		
		out_file = open('genes.fasta', 'w')
		out_file.write('>'+str(i)+'\n')
		out_file.write(seq1+'\n')
		out_file.write('>'+str(j)+'\n')
		out_file.write(seq2+'\n')
		out_file.close()
		
		muscle_cline = MuscleCommandline(input="genes.fasta", gapopen = -float(d))
		stdout, stderr = muscle_cline()

		align = AlignIO.read(StringIO(stdout), "fasta")
		print(align)
		print()
		alignments = pairwise2.align.globalds(Seq(seq1), Seq(seq2), blos, -d, -d)
		print('\t'+alignments[0][0][:disp] + '...'*int(len(alignments[0][0]) > disp))
		print('\t'+alignments[0][1][:disp] + '...'*int(len(alignments[0][1]) > disp))
		print('\tscore = '+str(alignments[0][2]))
	
		input(' ... ')
		










