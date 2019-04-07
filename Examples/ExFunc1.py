
from CMS import optimize as opt

def func1(x):
	"""
	func1 returns the distance between a
	point on the plane z = m1*x + m2*y + b
	and the origin of the coordinate system.
	"""
	
	m1 = 2
	m2 = -1
	b = 0
	
	z = m1*x[0] + m2*x[1] + b
	
	return (x[0]**2 + x[1]**2 + z**2)**.5
	return abs((x[0] - 1)**2 - x[0])
	
#func, init_x, max_x=None, min_x=None, tail=25, tmax=10**3, depth=2, init_dx=None, eps=1e-8

init_x = [-.5, 1.5]


x,func_val = opt.CMS(func1, init_x)

print(x)
print(func_val)
