import math
import matplotlib.pyplot as plt
from numpy.random import normal
from scipy.integrate import odeint
from SALib.analyze import fast

#Screening Vars
stRNA_bindings = normal(140., scale=10., size=20)
rf1_bindings = normal(34.4, scale=5., size=20)

stRNA_unbindings = normal(60.23, scale=5., size=20)
rf1_unbindings = normal(0.19, scale=0.05, size=20)

initial_codons = [200., 400., 1000., 2000.]

#Vars 
initial_codon_stRNA = 0.
initial_codon_RF1 = 0.

def sys(y, t):
     codon = y[0]
     codon_stRNA = y[1]
     codon_RF1 = y[2]
     #the model equations 
     d_codon_dt = (stRNA_unbind * codon_stRNA + rf1_unbind * codon_RF1) - ((stRNA_bind + rf1_bind) * codon)
     d_codon_stRNA_dt = (stRNA_bind * codon) - (stRNA_unbind * codon_stRNA)
     d_codon_RF1_dt = (rf1_bind * codon) - (rf1_unbind * codon_RF1)
     return [d_codon_dt, d_codon_stRNA_dt, d_codon_RF1_dt]

def solve(init_cond, time):
  result = odeint(sys, init_cond, time)
  codon_res = result[:, 0]
  codon_stRNA_res = result[:,1]
  codon_RF1_res = result[:, 2]

  ratio = [(codon_res[i], codon_stRNA_res[i], codon_RF1_res[i]) for i in range(len(codon_res))]

  #print(ratio)

  print('Ratio at Equilibrium: ' + str(ratio[len(ratio) - 1]))
  prob = math.pow((ratio[len(ratio) - 1][1]/initial_codon), (float(initial_codon/200.)))
  print('Probability of read through for ' + str(initial_codon/200) + ' stop codons: ' + str(prob))

  return prob, ratio

def run(initial_codon, initial_codon_stRNA, initial_codon_RF1):
  init_codon = initial_codon
  init_codon_stRNA = initial_codon_stRNA
  init_codon_RF1 = initial_codon_RF1
  init_cond = [init_codon,
              init_codon_stRNA,
              init_codon_RF1]
  time = list(range((3600))) #per minute for 24hrs
  return solve(init_cond, time)

def plot(time, to_plot):
  plt.plot(time, list([tup[0] for tup in to_plot]), label='codon')
  plt.plot(time, list([tup[1] for tup in to_plot]), label='codon:stRNA')
  plt.plot(time, list([tup[2] for tup in to_plot]), label='codon:RF1')
  plt.ylabel('Codons per State')
  plt.legend()
  plt.xlabel('Time (seconds)')
  plt.savefig('model.png')

results = []
params =[]
probabilities = []
counter = 0
for stRNA_bind in stRNA_bindings:
  for rf1_bind in rf1_bindings:
    for stRNA_unbind in stRNA_unbindings:
      for rf1_unbind in rf1_unbindings:
          for initial_codon in initial_codons:
            print('Model ' + str(counter) + ':')
            print('Screening Parameters: stRNA bind - ' + str(stRNA_bind) + ', RF1 bind - ' + str(rf1_bind) + ', stRNA unbind - ' + str(stRNA_unbind) + ', RF1 unbind - ' + str(rf1_unbind) +', codons - ' + str(initial_codon))
            res = run(initial_codon, initial_codon_stRNA, initial_codon_RF1)
            results.append(res[1])
            params.append((stRNA_bind, rf1_bind, stRNA_unbind, rf1_unbind, initial_codon))
            probabilities.append(res[0])
            counter = counter + 1

#plt.plot([1, 2, 5, 10], probabilities)
#plt.ylabel('Probability')
#plt.yscale('log')
#plt.legend()
#plt.xlabel('Number of Stop Codons')
#plt.savefig('prob_per_stop.png')

#Sensitivity Analysis - FAST
problem = {'num_vars': 5, 
           'names': ['stRNA_binding', 
           'rf1_binding', 
           'stRNA_unbinding', 
           'rf1_unbinding', 
           'initial_codons'], 
           'bounds': [[stRNA_bindings.min(), stRNA_bindings.max()], [rf1_bindings.min(), rf1_bindings.max()], 
           [stRNA_unbindings.min(), stRNA_unbindings.max()], [rf1_unbindings.min(), rf1_unbindings.max()], 
           [200., 2000.]]}

print('FAST')
Si = fast.analyze(problem, probabilities, print_to_console=False)
print('First order indices: ' + str(Si['S1']))
print('Total order indices: ' + str(Si['ST']))
