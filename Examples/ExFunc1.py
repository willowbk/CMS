
from CMS import optimize as opt

def func1(x):
	"""
	func1 returns the distance between
	a parabola and a line.
	"""
	
	return abs((x[0] - 1)**2 + 1 - x[0])
	
init_xs = [[-1 + .1*i] for i in range(20 + 1)]

for init_x in init_xs:
	

	x,func_val = opt.CMS(func1, init_x)

	print(init_x, x, func_vale)
	
