import itertools

with open('records/model_rec.txt') as f:
	results = f.readlines()
recs = []
for i in range(int(len(results)/6)):
	recs.append(results[i*6:6 + i*6])
recs = [(float(j[1]), float(j[2]), float(j[3]), float(j[4]), float(j[5][:-4])) for j in recs if float(j[5][:-4]) != 0.0]
grouped = list(set(sorted(recs, key=lambda k: (k[0], k[1], k[2], k[3]))))

imm_st = set([i[0] for i in grouped])
col_st = set([i[1] for i in grouped])
imm_rbs = set([i[2] for i in grouped])
col_rbs = set([i[3] for i in grouped])

imm_str_missing = [list(col_st), list(imm_rbs), list(col_rbs)]
col_str_missing = [list(imm_st), list(imm_rbs), list(col_rbs)]
imm_rbs_missing = [list(imm_st), list(col_st), list(col_rbs)]
col_rbs_missing = [list(imm_st), list(col_st), list(imm_rbs)]

imm_str_missing_perm = list(itertools.product(*imm_str_missing))
col_str_missing_perm = list(itertools.product(*col_str_missing))
imm_rbs_missing_perm = list(itertools.product(*imm_rbs_missing))
col_rbs_missing_perm = list(itertools.product(*col_rbs_missing))

imm_str_sensitivities = []
temp_g = []
for group in imm_str_missing_perm:
	temp_g = [(r[0], r[4]) for r in grouped if r[1] == group[0] and r[2] == group[1] and r[3] == group[2]]
	imm_str_sensitivities.append(temp_g)
	temp_g = []

col_str_sensitivities = []
for group in col_str_missing_perm:
	temp_g = [(r[1], r[4]) for r in grouped if r[0] == group[0] and r[2] == group[1] and r[3] == group[2]]
	col_str_sensitivities.append(temp_g)
	temp_g = []

imm_rbs_sensitivities = []
for group in imm_rbs_missing_perm:
	temp_g = [(r[2], r[4]) for r in grouped if r[0] == group[0] and r[1] == group[1] and r[3] == group[2]]
	imm_rbs_sensitivities.append(temp_g)
	temp_g = []

col_rbs_sensitivities = []
for group in col_rbs_missing_perm:
	temp_g = [(r[3], r[4]) for r in grouped if r[0] == group[0] and r[1] == group[1] and r[2] == group[2]]
	col_rbs_sensitivities.append(temp_g)
	temp_g = []

time_diff_imm = []
for compare in imm_str_sensitivities:
	comparisons = list(itertools.combinations(compare, 2))
	compare_group = []
	compare_gr_avg = []
	for comp in comparisons:
		str_diff = abs(comp[0][0] - comp[1][0])
		time_diff = abs(comp[0][1] - comp[1][1])
		diff_per_percent_str = time_diff/(str_diff * 100)
		compare_group.append(diff_per_percent_str)
	compare_gr_avg = sum(compare_group)/float(len(compare_group))
	time_diff_imm.append(compare_gr_avg)

avg_imm_st_change = sum(time_diff_imm)/len(time_diff_imm)
print('Average change to killswitch activation time per 1% of change in imm promoter strength: ' + str(avg_imm_st_change) + ' hrs') 

time_diff_col = []
for compare in col_str_sensitivities:
	comparisons = list(itertools.combinations(compare, 2))
	compare_group = []
	compare_gr_avg = []
	for comp in comparisons:
		str_diff = abs(comp[0][0] - comp[1][0])
		time_diff = abs(comp[0][1] - comp[1][1])
		diff_per_percent_str = time_diff/(str_diff * 100)
		compare_group.append(diff_per_percent_str)
	compare_gr_avg = sum(compare_group)/float(len(compare_group))
	time_diff_col.append(compare_gr_avg)

avg_col_st_change = sum(time_diff_col)/len(time_diff_imm)
print('Average change to killswitch activation time per 1% of change in col promoter strength: ' + str(avg_col_st_change) + ' hrs')

time_diff_rbs_imm = []
for compare in imm_rbs_sensitivities:
	if len(compare) == 1 or len(compare) == 0:
		continue
	comparisons = list(itertools.combinations(compare, 2))
	compare_group = []
	compare_gr_avg = []
	for comp in comparisons:
		str_diff = abs(comp[0][0] - comp[1][0])
		time_diff = abs(comp[0][1] - comp[1][1])
		diff_per_percent_aff = time_diff/(str_diff * 100)
		compare_group.append(diff_per_percent_aff)
	compare_gr_avg = sum(compare_group)/float(len(compare_group))
	time_diff_rbs_imm.append(compare_gr_avg)

avg_imm_rbs_change = sum(time_diff_rbs_imm)/len(time_diff_rbs_imm)
print('Average change to killswitch activation time per 1% of change in imm rbs affinity: ' + str(avg_imm_rbs_change) + ' hrs')

time_diff_rbs_col = []
for compare in col_rbs_sensitivities:
	if len(compare) == 1 or len(compare) == 0:
		continue
	comparisons = list(itertools.combinations(compare, 2))
	compare_group = []
	compare_gr_avg = []
	for comp in comparisons:
		str_diff = abs(comp[0][0] - comp[1][0])
		time_diff = abs(comp[0][1] - comp[1][1])
		diff_per_percent_aff = time_diff/(str_diff * 100)
		compare_group.append(diff_per_percent_aff)
	compare_gr_avg = sum(compare_group)/float(len(compare_group))
	time_diff_rbs_col.append(compare_gr_avg)

avg_col_rbs_change = sum(time_diff_rbs_col)/len(time_diff_rbs_imm)
print('Average change to killswitch activation time per 1% of change in col rbs affinity: ' + str(avg_col_rbs_change) + ' hrs')
