import systems
import rates
import itertools
import pickle
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.integrate import odeint

#param pairs to screen
copy_num_pairs = list(itertools.combinations(range(len(rates.copy_num)), 2)) + [tup[::-1] for tup in list(itertools.combinations(range(len(rates.copy_num)), 2))]+ list(zip(range(len(rates.copy_num)), range(len(rates.copy_num))))

promoter_pairs = list(itertools.combinations(range(len(rates.anderson_str)), 2)) +[tup[::-1] for tup in list(itertools.combinations(range(len(rates.anderson_str)), 2))] + list(zip(range(len(rates.anderson_str)), range(len(rates.anderson_str))))

rbs_pairs = list(itertools.combinations(range(len(rates.rbs_affinity)), 2)) + [tup[::-1] for tup in list(itertools.combinations(range(len(rates.rbs_affinity)), 2))] + list(zip(range(len(rates.rbs_affinity)), range(len(rates.rbs_affinity))))

#init imm only
init_imm_mRNA_sys0 = 0.
init_imm_sys0 = 0.
init_cond_sys0 = [init_imm_mRNA_sys0,
                  init_imm_sys0]
time_sys0 = list(range((24*60) - 1)) #24hrs

#init imm and col
init_col_mRNA_sys1 = 0.
init_col_sys1 = 0.
time_sys1 = list(range(24*60)) #24hrs

#init imm and col maxicell
time_sys2 = list(range((240*60) + 1)) #72hrs

#model count
num_models = 0

#model settings
settings = []

#run screen
for i in copy_num_pairs:
  for j in promoter_pairs:
    for k in rbs_pairs:
      model = systems.Systems(i[0], i[1], j[0], j[1], k[0], k[1])

      #solve imm only
      result_sys0 = odeint(model.imm_only, init_cond_sys0, time_sys0)
      imm_mRNA_sys0 = result_sys0[:, 0]
      imm_sys0 = result_sys0[:, 1]

      #init imm and col
      init_imm_mRNA_sys1 = imm_mRNA_sys0[len(imm_mRNA_sys0) - 1]
      init_imm_sys1 = imm_sys0[len(imm_sys0) - 1]
      init_cond_sys1 = [init_imm_mRNA_sys1,
                        init_imm_sys1,
                        init_col_mRNA_sys1,
                        init_col_sys1]

      #solve imm and col
      result_sys1 = odeint(model.imm_and_col, init_cond_sys1, time_sys1)
      imm_mRNA_sys1 = result_sys1[:, 0]
      imm_sys1 = result_sys1[:, 1]
      col_mRNA_sys1 = result_sys1[:, 2]
      col_sys1 = result_sys1[:, 3]

      #init imm and col maxicell
      init_imm_mRNA_sys2 = imm_mRNA_sys1[len(imm_mRNA_sys1) - 1]
      init_imm_sys2 = imm_sys1[len(imm_sys1) - 1]
      init_col_mRNA_sys2 = col_mRNA_sys1[len(col_mRNA_sys1) - 1]
      init_col_sys2 = col_sys1[len(col_mRNA_sys1) - 1]
      init_cond_sys2 = [init_imm_mRNA_sys2,
                        init_imm_sys2,
                        init_col_mRNA_sys2,
                        init_col_sys2]

      #solve imm and col maxicell
      result_sys2 = odeint(model.imm_and_col_maxicell, init_cond_sys2, time_sys2)
      imm_mRNA_sys2 = result_sys2[:, 0]
      imm_sys2 = result_sys2[:, 1]
      col_mRNA_sys2 = result_sys2[:, 2]
      col_sys2 = result_sys2[:, 3]

      #stitch together results of different systems
      imm_total = list(imm_sys0) + list(imm_sys1) + list(imm_sys2)
      col_total = [0.0 for i in range((24*60) - 1)] + list(col_sys1) + list(col_sys2)

      imm_mRNA_total = list(imm_mRNA_sys0) + list(imm_mRNA_sys1) + list(imm_mRNA_sys2)
      col_mRNA_total = [0.0 for i in range((24*60) - 1)] + list(col_mRNA_sys1) + list(col_mRNA_sys2)

      result_sys2_l = list(result_sys2)

      crossover = 0

      #maxicell active timeframe
      for l in range(len(result_sys2_l)):
        if result_sys2_l[l][3] > result_sys2_l[l][1]:
          crossover = l
          break

      print('Model: ' + str(num_models))
      print('Imm Copy Number: ' + str(rates.copy_num[i[0]]))
      print('Col Copy Number: ' + str(rates.copy_num[i[1]]))
      print('Imm Promoter Strength: ' + str(rates.anderson_str[j[0]]))
      print('Col Promoter Strength: ' + str(rates.anderson_str[j[1]]))
      print('Imm RBS Affinity: ' + str(rates.rbs_affinity[k[0]]))
      print('Col RBS Affinity: ' + str(rates.rbs_affinity[k[1]]))
      print('Maxicell Active Timeframe: ' + str(float(crossover)/60.) + ' hrs')
      print('------------------------------------------------------------------')

      settings.append((i[0], i[1], j[0], j[1], k[0], k[1]))

      num_models = num_models + 1

      file = open('records/model' + str(num_models-1) + 'imm_mRNA_total', 'wb')
      pickle.dump(imm_mRNA_total, file)
      file.close()
      file = open('records/model' + str(num_models-1) + 'col_mRNA_total', 'wb')
      pickle.dump(col_mRNA_total, file)
      file.close()
      file = open('records/model' + str(num_models-1) + 'imm_total', 'wb')
      pickle.dump(imm_total, file)
      file.close()
      file = open('records/model' + str(num_models-1) + 'col_total', 'wb')
      pickle.dump(col_total, file)
      file.close()

      file = open("records/model_rec.txt", 'a')
      file.write(str(num_models))
      file.write('\n')
      file.write(str(rates.anderson_str[j[0]]))
      file.write('\n')
      file.write(str(rates.anderson_str[j[1]]))
      file.write('\n')
      file.write(str(rates.rbs_affinity[k[0]]))
      file.write('\n')
      file.write(str(rates.rbs_affinity[k[1]]))
      file.write('\n')
      file.write(str(float(crossover)/60.) + ' hrs')
      file.write('\n')
      file.close()

#full timespan
time_total = list(range(288*60)) #288 hrs

#plot
fig = plt.figure(figsize=(10, 10))

#mRNA
plt.subplot(311)
for switch_sys in range(num_models):
  file = open('records/model' + str(switch_sys) + 'imm_mRNA_total', 'rb')
  to_plt = pickle.load(file)
  file.close()
  plt.plot(time_total, to_plt, label='imm mRNA, cn:' + str(settings[switch_sys][0]) + ', str:' + str(settings[switch_sys][2]) + ', rbs:' + str(settings[switch_sys][4]))
  file = open('records/model' + str(switch_sys) + 'col_mRNA_total', 'rb')
  to_plt = pickle.load(file)
  file.close()
  plt.plot(time_total, to_plt, label='col mRNA, cn:' + str(settings[switch_sys][1]) + ', str:' + str(settings[switch_sys][3]) + ', rbs:' + str(settings[switch_sys][5]))
plt.ylabel('mRNA Present')
plt.yscale('log')
xcoords = [(24*60) - 1 , 48*60]
for xc in xcoords:
    plt.axvline(x=xc, color = 'r', ls='dashed')

plt.text((24*60) - 1 + 0.1,  100, 'col transformation', rotation=90, color = 'r')
plt.text(48*60 + 0.1, 100, 'maxicell induction', rotation=90, color = 'r')

#protein
plt.subplot(312)
for switch_sys in range(num_models):
  file = open('records/model' + str(switch_sys) + 'imm_total', 'rb')
  to_plt = pickle.load(file)
  file.close()
  plt.plot(time_total, to_plt, label='imm, cn:' + str(settings[switch_sys][0]) + ', str:' + str(settings[switch_sys][2]) + ', rbs:' + str(settings[switch_sys][4]))
  file = open('records/model' + str(switch_sys) + 'col_total', 'rb')
  to_plt = pickle.load(file)
  file.close()
  plt.plot(time_total, to_plt, label='col, cn:' + str(settings[switch_sys][1]) + ', str:' + str(settings[switch_sys][3]) + ', rbs:' + str(settings[switch_sys][5]))
plt.ylabel('Protein Present')
plt.yscale('log')
for xc in xcoords:
    plt.axvline(x=xc, color = 'r', ls='dashed')

plt.savefig('model.png')
