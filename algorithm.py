import system as system
from time import time
import numpy as np

system.Import_Data('e_coli_ppm0_5-fold')
system.Set_Data('all')
system.Set_Params(19)

start_i = 0
t = 0
sol_fnd = False
tmax = 10**6
m = 25
eps = .05

seq = system.Generate_Seq()
F = system.Get_Fitness(seq)

best_seq = list(seq)
Fs = []
Ff = 0

path = '_'.join(seq)
init_seq = list(seq)

best_F = F
start_F = F
print('START: ' + str(best_F) + ' ... ' + str(seq))
F_t = []
seq_t = []

system.dx = 10
system.depth = 2
start_time = time()

while t < tmax:
	
	nn_ws = system.Get_All_Neighbors(seq)
	max_F = F
	nn_ws_count = 0
	i = start_i
	prev_start_i = start_i
	
	while max_F == F and nn_ws_count < len(nn_ws) and (Ff >= 0 or len(Fs) < m):
		nn_w = nn_ws[i]
		max_F = max([system.Get_Fitness(nn_w[0]), max_F])
		nn_ws_count += 1
		i = i + 1
		if i >= len(nn_ws):
			i = 0
		t += 1
		
		if len(Fs) == m:
			Fs = Fs[1:m + 1]
		Fs.append(F)
		
		Ff = (Fs[-1] - Fs[0])*(tmax - t)/len(Fs) + F 
	
	start_i = i - 1 + int(nn_ws_count == len(nn_ws))
	max_nn = list(nn_w[0])
	
	if max_F <= F:
		
		system.dx /= 2
		if system.dx < 1.e-4:
			system.dx = 10
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
			new_F = system.Get_Fitness(new_seq)
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
	
	if F > best_F:
		print(t, 'BEST: ' + str(system.Get_Fitness(seq)) + ' ... ' + str(seq))
		
		best_seq = list(seq)
		best_F = F
		
		F_t.append([F, t])
		seq_t.append([float(term) for term in seq] + [t])
		
runtime = time()-start_time
print("Completion time: " + str(round(runtime, 3)) + " seconds")
print('dx = '+str(system.dx))
print('depth = '+str(system.depth))
print('BEST: ' + str(best_F) + ' ... ' + str(best_seq))


with open("CodonBias_Results_"+organism+"_"+str(data_set)+".out", "a") as out_file:
	out_file.write('START: ' + str(start_F) + ' SEQi: ' + str(init_seq) + ' BEST: ' + str(best_F) + ' SEQf: ' + str(best_seq) + ' TIME: '  + str(round(runtime, 3)) + ' seconds fin. dep. = ' + str(system.depth) + '\n')
	
	
	
