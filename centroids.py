'''
Plotting CA, mutation map and mutational age.

For multiple biopsies plot allele frequency and detected alleles

detection_threshold tunable, as is biopsy_num and biopsy radius (r)

Jacob Scott 23 December 2015, updated 21 Jan 2016 with new algorithm for biopsy

usage: python centroids.py read_path write_path filename size time 

'''
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from pylab import rcParams
import os
import clonal_evolution_functions as cef
from scipy.spatial import distance
import scipy.stats as stats
import sys

# Parse input arguments
read_path = str(sys.argv[1]) 
write_path = str(sys.argv[2]) 
filename = str(sys.argv[3])
size = int(sys.argv[4])
time = int(sys.argv[5])

SI1 = 0 #placeholders for Shannon Index values
colony_size_threshold = 1000 #what size colonies to consider in size distribution plot

my_cmap = plt.cm.get_cmap('nipy_spectral')
my_cmap.set_under('w')

#mutation flag data
data = open(read_path+'carriedMutation'+str(time)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('int')
CM1 = np.reshape(CA, (size,size))
N1 = np.count_nonzero(CA)
for species in range (1, np.amax(CA)): SI1 = SI1 + cef.shannon(np.bincount(CA)[species],N1)
SItrunc = float("{0:.4f}".format(SI1))

#population data
data = open(read_path+filename+'.txt').read().replace(',',' ').replace('\n',' ')
x = data.split()
cellmatrix = np.array(x).astype(str)
y = len(cellmatrix)/5
data_shaped = np.reshape(cellmatrix, (y,5))
cells_over_time = data_shaped[:,2:].astype(int) #grabs stem, non-stem, and associated timesteps
muts_over_time = data_shaped[:,1].astype(int) #grabs unique mutation events (children)
total_over_time = np.zeros(len(cells_over_time))
for i in range(0,len(cells_over_time)):
	total_over_time[i] = (cells_over_time[i,0] + cells_over_time[i,1])


RMSvalues = np.zeros(np.max(CM1))
distances_to_centroid = []
colony_size_list = []

# calculate distance to centroids of each colony and count colony sizes
for mutation in range(1, np.max(CM1)):
	x,y = np.where(CM1 == mutation)
	centroid = (np.mean(x), np.mean(y))
	muts = list(zip(x,y))
	colony_size_list.append(len(muts))
	delta_list = []
	total_distance = 0
	for point in muts: 
		total_distance += distance.euclidean(point,centroid)**2
		RMSvalues[mutation] = total_distance/len(muts)
		delta_list.append(distance.euclidean(point,centroid))
	distances_to_centroid.append(delta_list)

#delete columns with colonies less than colony_size_threshold cells
distances_to_centroid_trunc = []
trunc_list = np.zeros(len(distances_to_centroid), dtype = bool)
i = 0
for thing in distances_to_centroid:
	if len(thing)>colony_size_threshold:
		trunc_list[i] = 1
	i += 1

CSD_mean = np.mean(colony_size_list)
CSD_std = np.std(colony_size_list)
CSD_skew = stats.skew(colony_size_list)
stds = np.zeros(len(distances_to_centroid))
skews = np.zeros(len(distances_to_centroid))
mut_labels = np.linspace(0, len(distances_to_centroid), len(distances_to_centroid)).astype('int')

#find moments of distribution
for i in range(0,len(distances_to_centroid)):
	stds[i] = np.std(distances_to_centroid[i])
	skews[i] = stats.skew(distances_to_centroid[i])

mean_distances = np.zeros(len(distances_to_centroid))
for i in range(0,len(distances_to_centroid)): mean_distances[i] = np.mean(distances_to_centroid[i])
info_array = np.vstack((mean_distances, distances_to_centroid, stds, skews, mut_labels))

info_array = info_array[:,trunc_list==1] #slice array so that only columns in trunc_list are kept

# sort by descending means
idx = np.argsort(info_array[0])
info_array=info_array[:,idx]

# sort together 
sorted_distances = []
for element in info_array[1]: sorted_distances.append(element)
sorted_distances.reverse()
stds = info_array[2]
skews = info_array[3]
mut_labels = info_array[4]
rev_skews = skews[::-1]
new_skews = np.zeros(len(stds)+1) 
rev_stds = stds[::-1]
new_stds = np.zeros(len(stds)+1) 
for i in range (0,len(stds)): 
	new_stds[i+1]=rev_stds[i]
	new_skews[i+1]=rev_skews[i]

'''PLOT BW of centroid_dist and STDEV and SKEW'''
rcParams['figure.figsize'] = 20,10
# plt.figure()
fig, ax1 = plt.subplots()
ax1.boxplot(sorted_distances)
ax1.set_xlabel('Mutation')
ax1.set_xticklabels(mut_labels)
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('Distance to centroid', color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')


ax2 = ax1.twinx()
ax2.plot(new_stds, 'r--', label = 'stds')
ax2.plot(new_skews, 'k--', label='skews')
ax2.set_ylabel('moments of variation', color='r')
ax2.legend()
for tl in ax2.get_yticklabels():
    tl.set_color('r')
plt.savefig(write_path+filename+'_BW_centroid_dist.png', dpi = 500)

# ''' PLOT POPULATION DYNAMICS'''
# plt.figure()
# plt.plot(cells_over_time[:,2],cells_over_time[:,1], color = 'k', linewidth = 5, label = 'non-stem')
# plt.plot(cells_over_time[:,2],cells_over_time[:,0], color = 'b', linewidth = 5, label = 'stem')
# plt.plot(cells_over_time[:,2],total_over_time, color = 'r', linewidth = 5, label = 'total')
# plt.plot(cells_over_time[:,2],muts_over_time, color = 'c', linewidth = 5, label = 'mutations')
# plt.yscale('log')
# plt.ylabel('Cells')
# plt.xlabel('Time')
# plt.legend(loc = 'upper left')
# plt.savefig(write_path+filename+'_pop_mut_dynamics.png', dpi = 500)


# ''' PLOT COLONY SIZE DISTRIBUTION '''
# plt.figure()
# plt.hist(colony_size_list, bins = np.logspace(0, 6), normed = 0)
# plt.gca().set_xscale('log')
# plt.title('Colony size distribution. Mean: '+str(CSD_mean)+' stdev: '+str(CSD_std)+' skew: '+str(CSD_skew))
# plt.xlabel('Colony size')
# plt.ylabel('Frequency')
# plt.savefig(write_path+filename+'colony_size_dist'+'.png', dpi = 500)

# ''' PLOT CA '''
# rcParams['figure.figsize'] = 10,7
# plt.figure()
# plt.pcolor(CM1, cmap='nipy_spectral', vmin = 0.001)
# plt.title('CA_Shannon:'+str(-SItrunc))
# plt.colorbar()
# plt.savefig(write_path+filename+'_CA.png', dpi = 500)

''' PLOT distributions of moments of variation across whole simulation '''
rcParams['figure.figsize'] = 10,20

skew_list = info_array[3].astype(float)
np.savetxt(write_path+'skews.txt', skew_list)
std_list = info_array[2].astype(float)
np.savetxt(write_path+'stds.txt', std_list)
mean_list = info_array[0].astype(float)
np.savetxt(write_path+'means.txt', mean_list)

moment_lists = [skew_list, std_list, mean_list]
data_labels = ['skew', 'std', 'mean']
plt.subplot(3,1,1)
plt.boxplot(moment_lists[0])
# plt.xticks([1], data_labels[0])
plt.title('Skews')

plt.subplot(3,1,2)
plt.boxplot(moment_lists[1])
# plt.xticks([1], data_labels[1])
plt.title('StdDevs')

plt.subplot(3,1,3)
plt.boxplot(moment_lists[2])
# plt.xticks([1], data_labels[2])
plt.title('Means')

plt.savefig(write_path+filename+'_momentDist.png', dpi = 500)

# plt.show()
