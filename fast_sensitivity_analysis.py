from SALib.analyze import fast
import numpy as np

with open('records/model_rec.txt') as f:
	results = f.readlines()
recs = []
for i in range(int(len(results)/6)):
	recs.append(results[i*6:6+i*6])
recs = [(float(j[1]), float(j[2]), float(j[3]), float(j[4]), float(j[5][:-4])) for j in recs if float(j[5][:-4]) != 0.0]

Y = np.array([k[4] for k in recs])
param_values = np.array([(l[0], l[1], l[2], l[3]) for l in recs])
problem = {'num_vars': 4, 'names': ['imm_st', 'col_st', 'imm_rbs', 'col_rbs'], 'bounds': [[0.01, 1.], [0.01, 1.], [0.07, 0.6], [0.07, 0.6]]}

print('FAST')
Si = fast.analyze(problem, Y, print_to_console=False)
print('First order indices: ' + str(Si['S1']))
print('Total order indices: ' + str(Si['ST']))
