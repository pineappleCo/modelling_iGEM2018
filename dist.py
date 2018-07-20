import matplotlib.pyplot as plt

with open('records/model_rec.txt') as f:
	recs = f.readlines()

times = []
for i in range(len(recs) - 1)[1:]:
	try:
	        times.append(recs[(i*6) - 1])
	except(IndexError):
		break

num_times = [float(t[:-6]) for t in times]
print(num_times.index(max(num_times)))
print(num_times.index(min([a for a in num_times if a != 0.])))

print(max(num_times))
print(num_times[456])

plt.hist(num_times)
plt.ylabel('count')
plt.xlabel('killswitch activation time (hrs)')
plt.savefig('model_dist.png')
