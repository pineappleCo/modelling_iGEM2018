from scipy.integrate import odeint

#variables - per 60min generation
rf1_mut = 1.99e-10
trna_mut = 1.99e-10
rf1_sup = 1.99e-10
trna_sup = 1.99e-10
read_through_rt = 1.

def sys(y, t):
	normal_cells = y[0]
	mutant_cells = y[1]
	#the model equations
	d_normal_dt = ((rf1_sup * trna_sup) * mutant_cells) - ((rf1_mut * trna_mut) * normal_cells)
	d_mutant_dt = ((rf1_mut * trna_mut) * normal_cells) - ((rf1_sup * trna_sup) * mutant_cells)
	d_read_through_dt = read_through_rt * mutant_cells
	return [d_normal_dt, d_mutant_dt, d_read_through_dt]

#init
init_normal_cells = 10e14
init_mutant_cells = 0.
init_read_through = 0.
init_cond = [init_normal_cells, init_mutant_cells, init_read_through]
time = list(range(24))

#solve
result = odeint(sys, init_cond, time)
read_through_sys = result[:, 2]

#results
fail_rt = (read_through_sys[len(read_through_sys) - 1]/ 10e14) * 100
print("Failure Rate: " + str(fail_rt))
