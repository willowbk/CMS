
def Single_Moves(moves):
	"""
	    This function generates a set of all
	possible unique moves in which a single
	parameter is changed in either the 
	positive or negative direction.
	
	Parameters
	----------
		moves : list of chars
		        1-D array of moves represented as characters taking only integer values.
	Returns
	-------
		nns : set of strings
		        A set of all single-move changes to the move list.
	
	"""
	nns = set([])

	for i in range(len(moves)):
		nns.add('_'.join(moves[:i] + [str(int(moves[i]) + 1)] + moves[i+1:]))
		nns.add('_'.join(moves[:i] + [str(int(moves[i]) - 1)] + moves[i+1:]))
		
	return nns



def Get_All_Neighbors(x, dx, depth, min_x, max_x):
	
	"""
	    Given that single moves can be generated,
	this function will then construct a set
	of all combinations of single moves 
	up to the value of the parameter 'depth.'
	
	    Once the full set of moves are found,
	these moves are applied to the current 
	location 'x' using the step size 'dx'
	while remaining within the boundaries 
	specified by 'min_x' and 'max_x.'
	
	Parameters
	----------
		x : list of floats
			An array of floats representing the current location
			in parameter space.
		
		dx : float
			The step size for each move.

		depth : int
			How many single moves to apply for updates to 'x.'
		
		min_x : list of floats
			Lower bound on values of 'x.'
		
		max_x : list of floats
			Upper bound on values of 'x.'
		
		        
	Returns
	-------
		neighbors : list
		        A list of all possible modifications of 'x' with step size 'dx' and 
		        number of changes 'depth.'
	
	"""
	
	ns = Single_Moves(['0' for i in range(len(x))])

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
	
	neighbors = [[min([max([float(x[i]) + dx*int(move), min_x[i]]), max_x[i]]) for i,move in enumerate(n.split('_'))] for n in ns]
	
	return neighbors



def CMS(func, init_x, max_x=None, min_x=None, tail=25, tmax=10**3, depths=None, init_dx=None, eps=1e-8):
	
	"""
	    This algorithm minimizes the function 'func' on 
	the parameter space described by 'x' starting from
	'init_x.' To optimize, this algorithm attempts changing the
	current values of 'x' in all combinations described by the 
	'depths' list (by single changes, double changes, etc.) The 
	magnitude of each change is 'dx' regardless of the entry of 
	'x' being updated:
	
	Example:
		Given that 
			x1, x2, ...
		
		are the entries of 'x,' all single changes would be
		
			x1 + dx, x2, ...
			x1 - dx, x2, ...
			x1,      x2 + dx, ...
			etc.
		
		This set of moves is one the default settings 
		specified by the first entry of 'depths' being
		set to 1.
		
		All double changes (second entry of 'depths' set
		to 2 in the default settings) are given by
		
			x1 + 2dx, x2, ...
			x1 + dx, x2 + dx, ...
			x1 - dx, x2 + dx, ...
			etc.
	
	
	    After a full list of possible changes has been computed,
	the value of 'func' is evaluated for each change sequentially
	until an improvement in 'func' is found. If a certain move 
	further minimizes 'func,' this move is applied repeatedly until 
	no further improvement is found. To optimize this step, a binary
	search is employed to determine the appropriate number of 
	application of the minimizing move.
	
	    Throughout this algorithm, an average rate of convergence is 
	computed given the optimal value of 'func' several function 
	calls in the past (the number of calls is given by the parameter 
	'tail'). This average is then used to determine if the algorithm 
	is expected to find a minimum in 'func' given the remaining allotted
	function calls specified by 'tmax.' If not (or the full list of 
	moves is evaluated with an improvement found), the value of 'dx' is 
	updated, and the algorithm continues.
	
	
	Parameters
	----------
		func : callable func(x,*args)
			The objective function to be minimized.
		
		init_x : list of floats
			Initial values for each entry of 'x'

		max_x : int
			Upper bound on values of 'x.'
		
		min_x : list
			Lower bound on values of 'x.'
		
		tail : list
			Number of calls to 'func' before
			the average rate of convergence 
			is calculated.
		
		tmax : integer
			Total number of function calls
			to 'func' before termination.
		
		depths : list of integers
			Each entry specifies the number of 
			changes to be applied to the values of 
			'x' each of size 'dx.' 
			
		init_dx : float
			Initial value for the variable
			step size 'dx.'
		
		eps : float
			Minimum value of 'dx' below which
			the step size is reset to 'init_dx.'
			
		
		        
	Returns
	-------
		best_x : list
			The recovered list of values of 'x' which 
			most-minimized the objective function 'func.'
		        
		best_F : float
			Value of the objective function 'func'
			at the optimal parameter values 'best_x.'
	
	"""
	
	# Default ranges based roughly on Python precision.
	if max_x == None:
		max_x = [1e16 for i in range(len(init_x))]
		
	if min_x == None:
		min_x = [-1e16 for i in range(len(init_x))]
		
	if init_dx == None:
		# If not initial value of 'dx' is specified,
		# this value will default to the distance to
		# the nearest boundary.
		init_dx = min([min([max_x[i] - init_x[i], init_x[i] - min_x[i]]) for i in range(len(init_x))])
	
	if depths == None:
		# 'depths' defaults to 1 and 2 changes
		# unless there is only one value in 'x'
		# to be updated.
	
		if len(init_x) == 1:
			depths = [1]
		else:
			depths = [1,2]
	
	dx = init_dx
	start_i = 0
	t = 0
	
	x = list(init_x)
	F = func(x)
	Ff = 0
	Fs = []

	best_x = list(x)
	best_F = F
	start_F = F
	#print('START: ' + str(best_F) + ' ... ' + str(x))

	while t < tmax:
		
		nns = []
		for depth in depths:
			nns = nns + Get_All_Neighbors(x, dx, depth, min_x, max_x)
		
		min_F = F
		nns_count = 0
		i = start_i
		prev_start_i = start_i
		
		
		while min_F == F and nns_count < len(nns) and (Ff <= 0 or len(Fs) < tail):
			nn = nns[i]
			min_F = min([func(nn), min_F])
			nns_count += 1
			i = i + 1
			if i >= len(nns):
				i = 0
			t += 1
		
			if len(Fs) == tail:
				Fs = Fs[1:tail + 1]
			Fs.append(F)
		
			Ff = (Fs[-1] - Fs[0])*(tmax - t)/len(Fs) + F 
	
		start_i = i - 1 + int(nns_count == len(nns))
		min_nn = list(nn)
		
		if min_F >= F:
		
			dx /= 2
			if dx < eps:
				dx = init_dx
			Ff = 0
			Fs = []
		else:
			move = [min_nn[i] - xi for i,xi in enumerate(x)]
			x = list(min_nn)
			F = min_F
		
			itera = int(min([min([max_x[i] - x[i], x[i] - min_x[i]]) for i in range(len(init_x))])/(max(depths)*dx))
		
			while itera >= 1:
				new_x = [min([max([x[i] + itera*move[i], min_x[i]]), max_x[i]]) for i in range(len(x))]
				new_F = func(new_x)
				t += 1
			
				if new_F < F:
					x = list(new_x)
					F = new_F
					print('\t F = '+str(F), itera)
				else:
					itera /= 2
				
				Fs.append(F)
	
		if F < best_F:
			#print(t, 'BEST: ' + str(func(x)) + ' ... ' + str(x))
		
			best_x = list(x)
			best_F = F
		
	return best_x, best_F




	
	
