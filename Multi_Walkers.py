import system as system
from time import time
import numpy as np
import sys
import Analyze_CodonBias_Params as ACP

organism = ACP.organism
#system.Import_Data(organism + '_ppm0_fake_5-fold_19_param')
#system.Import_Data(organism + '_ppm0_fake_5-fold_16_param')
#system.Import_Data(organism + '_ppm0_fake_5-fold_12_param')
#system.Import_Data(organism + '_ppm0_fake_5-fold_7_param')
#system.Import_Data(organism + '_ppm0_fake_5-fold_3_param')

if ACP.data_params == 0:
	system.Import_Data(organism + '_ppm0_5-fold')
else:
	system.Import_Data(organism + '_ppm0_fake_5-fold_'+str(ACP.data_params)+'_param')

if len(sys.argv[1]) == 1:
	data_set = int(sys.argv[1])
else:
	data_set = sys.argv[1]
	
system.Set_Data(data_set)
system.Set_Params(ACP.model_params)

try:
	import Codon_Seqs as CS
	seq = CS.seqs[data_set]
	
except ImportError as err:
	print(err)
	seq = system.Generate_Seq()
#if str(data_set) == '0':
#	seq = ['5.5877073844080885', '3.5506165681886435', '12.013009058212262', '6.261984191869559', '4.047156681068582', '18.368079991895513', '0.0', '37.63786979538432', '0.000152587890625', '0.008402233274148237', '57.668928538555306', '0.028942105948472234', '0.13623338315301178', '60.923605431346104', '1.5542433154575566', '91.51531125148357', '61.605755716855874', '0.05758942773713036', '78.65886418833766', '0.0']
#
#elif str(data_set) == '1':
#	seq = ['5.346124019933182', '3.4535829117211634', '10.528006396245788', '5.653557420786928', '5.903932311801209', '17.463220924896035', '0.068359375', '43.651422269627105', '0.0', '0.0', '53.52080782204495', '0.0', '0.11590563293758893', '54.80318786232256', '1.464782440006682', '84.11187559946501', '54.83819407334332', '0.04349816118724818', '68.91798243187756', '0.0']
#
#elif str(data_set) == '2':
#	seq = ['5.540443161137184', '3.532324154144817', '10.636334365147505', '6.195928954537513', '5.276031818222799', '19.409667025694898', '0.0', '42.666731440411795', '0.0', '0.00590471223822597', '58.456528213810756', '0.03057560676146647', '0.1884417501542057', '60.4185668409187', '1.644416009927687', '91.52912185045383', '61.617062754143', '0.06588394090473601', '73.51823098181785', '0.0011678948213819332']
#
#elif str(data_set) == '3':
#	seq = ['5.573092287725981', '3.618686220961282', '10.829718098695453', '6.230678093680992', '5.101180003370167', '19.602359243961015', '0.0', '46.570532234365814', '0.0', '0.0', '61.64722377201966', '0.0017189073887274435', '0.15582890979262617', '61.93213364768916', '1.6878954629804208', '86.26026904680145', '60.19129668881764', '0.05354398809322926', '76.79628648400673', '0.0']
#
#elif str(data_set) == '4':
#	seq = ['5.010371912184767', '3.041012541742228', '10.508303155038922', '5.426833862803007', '4.701617901115666', '16.493605671866142', '0.0', '37.60716006285758', '0.0018310546875', '0.0', '48.79064776790922', '0.00091552734375', '0.09113453679454266', '53.78401411801228', '1.2253410682626618', '78.15876262859254', '56.798629182903206', '0.05716484156041038', '64.10101304784256', '0.000762939453125']
#seq = [3.4793115162811032, 0.77204196708050021, 1.4184851677091148, 6.8083811400534984, 4.263736104128478, 2.0689098614003436, 10.006395184647692, 0.23559291303599536, 0.0056411000307789311, 3.2892647938308457, 0.0, 0.0, 7.2987985768297676, 0.15628451455868056, 99.905169614062132, 36.830510594615177, 1.6623364795972033, 3.6432940453429801, 18.658971991099193]
#seq = [3.5363338419758805, 0.91283793571168892, 3.5679087165940992, 6.8772790578431806, 4.1800450674245466, 2.4301205427123351, 4.9821559793404857, 0.085603472441813211, 0.0, 2.4835462657059564, 10.481044475846566, 0.0, 6.7832992557607614, 0.0770550620386765, 6581.6991817025391, 69.745521839478755, 2.7491526135487465, 3.9572315725681655, 2463.6387808954291]
#seq = [str(term) for term in seq]


score_i = ACP.score_i
start_i = 0
t = 0
sol_fnd = False
tmax = 1000*1000
m = 25

F = -1
eps = .05

seq = [str(float(term)*(1.00 - eps + 2*eps*np.random.random())) for term in seq]
F = system.Get_Fitness(seq)[score_i]

	


best_seq = list(seq)
Fs = []
Ff = 0

path = '_'.join(seq)
init_seq = list(seq)
#F = system.Get_Fitness(seq, 'training')[score_i]
best_F = F
start_F = F
print('START: ' + str(best_F) + ' ... ' + str(seq))
mode = 'FD'


F_t = []
seq_t = []

system.dx = 10
system.depth = 2
start_time = time()

while t < tmax:
	
	if mode == 'BD':
		print('Getting neighbors ...')
		nn_ws = system.Get_All_Neighbors(seq)
		Fs = []
		print('Evaluating fitnesses ...')
		print('dx = '+str(system.dx))
		print('depth = '+str(system.depth))
		for nn_w in nn_ws:
			Fs.append(system.Get_Fitness(nn_w[0], 'training')[score_i])
			t += 1
	
		max_F = max(Fs)
		max_nn_w = nn_ws[Fs.index(max_F)]
		
		print('Maximum found.')
		while '_'.join(max_nn_w[0]) == path:
			#print('Maximum found.')
			Fs.remove(max_F)
			nn_ws.remove(max_nn_w)
		
			max_F = max(Fs)
			max_nn_w = list(nn_ws[Fs.index(max_F)])
			
		max_nn = max_nn_w[0]
		
	elif mode == 'FD':
		#print('Getting neighbors ...')
		nn_ws = system.Get_All_Neighbors(seq)
		max_F = F
		nn_ws_count = 0
		#print('Evaluating fitnesses ...')
		#print('dx = '+str(system.dx))
		#print('depth = '+str(system.depth))
		#print(max_F)
		i = start_i
		prev_start_i = start_i
		
		while max_F == F and nn_ws_count < len(nn_ws) and (Ff >= 0 or len(Fs) < m):
			nn_w = nn_ws[i]
			max_F = max([system.Get_Fitness(nn_w[0])[score_i], max_F])
			nn_ws_count += 1
			i = i + 1
			if i >= len(nn_ws):
				i = 0
			t += 1
			
			if len(Fs) == m:
				Fs = Fs[1:m + 1]
			Fs.append(F)
			
			Ff = (Fs[-1] - Fs[0])*(tmax - t)/len(Fs) + F 
			#print(Ff, len(Fs))
		start_i = i - 1 + int(nn_ws_count == len(nn_ws))
		max_nn = list(nn_w[0])
		#print(str(prev_start_i) + ' + ' + str(nn_ws_count - 1) + ' = ' + str(start_i))
		
	elif mode == 'FS':
		max_F = F
		nn_ws_count = 0
		print('Evaluating fitnesses ...')
		
		while max_F == F and nn_ws_count < 100:
			nn = system.Update_Seq(seq)
			max_F = max([system.Get_Fitness(nn)[score_i], max_F])
			nn_ws_count += 1
			#print(' ... ' + str(max_F))
			print(nn_ws_count)
		max_nn = list(nn)
		
	
	#print('Comparing to current fitness ...')
	if max_F <= F:
		#print('Only deleterious or neutral moves. Changing grain.')
		system.dx /= 2
		if system.dx < 1.e-4:
			system.dx = 10
			#if not sol_fnd:
			#	system.depth += 1
			sol_fnd = False
		Ff = 0
		Fs = []
	else:
		move = [float(max_nn[i]) - float(seq[i]) for i in range(len(seq))]
		seq = list(max_nn)
		F = max_F
		
		itera = 10**6
		
		while itera >= 1:
			new_seq = [str(min([max([float(seq[i]) + itera*move[i], system.min_vals[i]]), system.max_vals[i]])) for i in range(len(seq))]
			new_F = system.Get_Fitness(new_seq)[score_i]
			t += 1
			
			if new_F > F:
				seq = list(new_seq)
				F = new_F
				print('\t F = '+str(F), itera)
			else:
				itera /= 10
				
			Fs.append(F)
		
		path = '_'.join(seq)
		sol_fnd = True
	
	#print('CURR: ' + str(F) + ' ... ' + str(seq))
	if F > best_F:
		print(t, 'BEST: ' + str(system.Get_Fitness(seq)) + ' ... ' + str(seq))
		#dseq = [float(seq[i]) - float(best_seq[i]) for i in range(len(seq))]
		#for j in range(1, 100 + 1):
		#	proj_seq = [str(max([float(best_seq[i]) + (2*j/100)*dseq[i], 0])) for i in range(len(seq))]
		#	proj_F = system.Get_Fitness(proj_seq, 'training')[score_i]
		#	t += 1
		#	if proj_F > F:
		#		F = proj_F
		#		seq = list(proj_seq)
		#		print('projected', F, j)
		
		#print('t = '+str(t))
		best_seq = list(seq)
		best_F = F
		
		F_t.append([F, t])
		seq_t.append([float(term) for term in seq] + [t])
		
runtime = time()-start_time
print("Completion time: " + str(round(runtime, 3)) + " seconds")
print('dx = '+str(system.dx))
print('depth = '+str(system.depth))
print('BEST: ' + str(best_F) + ' ... ' + str(best_seq))

if False:
	F_t = list(np.array(F_t).T)
	seq_t = list(np.array(seq_t).T)
	
	out_file = open('Multi_Walkers_codon.py', 'w')
	
	out_file.write('F_t = '+str([list(ar) for ar in F_t])+'\n')
	out_file.write('seq_t = '+str([list(ar) for ar in seq_t]))
	out_file.close()


	exit()

with open("CodonBias_Results_"+organism+"_"+str(data_set)+".out", "a") as out_file:
	out_file.write('START: ' + str(start_F) + ' SEQi: ' + str(init_seq) + ' BEST: ' + str(best_F) + ' SEQf: ' + str(best_seq) + ' TIME: '  + str(round(runtime, 3)) + ' seconds fin. dep. = ' + str(system.depth) + '\n')
	
	
	
