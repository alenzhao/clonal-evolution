'''
Plotting CA, mutation map and mutational age.

For multiple biopsies plot allele frequency and detected alleles

detection_threshold tunable, as is biopsy_num and biopsy radius (r)

Jacob Scott 23 December 2015, updated 21 Jan 2016 with new algorithm for biopsy

'''
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from pylab import rcParams
import os
import clonal_evolution_functions as cef
from scipy.spatial import distance

# parameters
size = 500 #size of the array
time = 500
SI1 = 0 #placeholders for Shannon Index values


read_path = '../andrea_test/non-stem/text/'
write_path = '../andrea_test/non-stem/figs/'

#mutation flag data
data = open(read_path+'carriedMutation'+str(time)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('int')
CM1 = np.reshape(CA, (size,size))
N1 = np.count_nonzero(CA)
for species in range (1, np.amax(CA)): SI1 = SI1 + cef.shannon(np.bincount(CA)[species],N1)
SItrunc = float("{0:.4f}".format(SI1))

centroid_dict = {}
RMSvalues = np.zeros(np.max(CM1))

for mutation in range(1, np.max(CM1)):
	x,y = np.where(CM1 == mutation)
	centroid = (np.mean(x), np.mean(y))
	muts = list(zip(x,y))
	centroid_dict.update({mutation: centroid})
	total_distance = 0
	for point in muts: total_distance += distance.euclidean(point,centroid)**2
	RMSvalues[mutation] = total_distance/len(muts)

'''PLOTS'''
rcParams['figure.figsize'] = 10,10
plt.figure()


plt.hist(RMSvalues, alpha=0.4)
# plt.xticks(allele_ticks, rotation = 315)
plt.ylabel('frequecy')	
plt.xlabel('RMS value')
# plt.ylim([0, 1])
plt.title('Histogram of RMS values')

plt.figure()
plt.pcolor(CM1)
plt.title('CA_Shannon_'+str(-SItrunc))
# plt.savefig(write_path+str(biopsy_num)+'TEST2 r= '+str(r)+'.png', dpi = 500)

plt.show()
