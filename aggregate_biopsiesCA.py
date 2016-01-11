'''
Plotting CA, mutation map and mutational age.

Jacob Scott 23 December 2015

'''
import matplotlib.pyplot as plt
import numpy as np
import random
from scipy.spatial import distance
from collections import Counter
from matplotlib.patches import Circle
from math import log as ln
from pylab import rcParams

#calculate the shannon index
def shannon(n, N):
    """ Relative abundance """
    if n == 0:
        return 0
    else:
        return (float(n)/N) * ln(float(n)/N)

#sum up the digits in a string
def sum_digits(digit):
    return sum(int(x) for x in digit if x.isdigit())

def gather_biopsies(biopsy_num, r):
	while len(biopsy_sites) < biopsy_num:
		# newpoint = [random.randint(r,size-r),random.randint(r,size-r)] #not including over the edge
		newpoint = [random.randint(0,size),random.randint(0,size)] #overlap OK, over edge ok
		biopsy_sites.append(newpoint)
	return biopsy_sites

def do_biopsies_aggregate(size, biopsy_num, r, CM1, biopsy_sites):
	area = 4*r**2
	biopsy_Mutlist = np.zeros((biopsy_num,area)).astype('int')
	aggregate_biopsy = np.array([]).astype('int')
	cell_count_inBx = np.zeros(biopsy_num)
	SIBx_agg = 0
	for row in range(0, size):
		for column in range(0, size):
			a = (row,column)
			for bx in range(0, biopsy_num):
				punch = biopsy_sites[bx]
				if distance.euclidean(a,punch) <= r:
					biopsy_Mutlist[bx][cell_count_inBx[bx]] = CM1[column][row]
					cell_count_inBx[bx] += 1 
					aggregate_biopsy = np.append(aggregate_biopsy,CM1[column][row])
	for bx in range(0, biopsy_num):
		SIBx_temp = 0
		biopsy_Mutlist_temp = (biopsy_Mutlist[bx])[0:cell_count_inBx[bx]]
		# print biopsy_Mutlist_temp
		for x in range (0, np.amax(biopsy_Mutlist_temp)):
			SIBx_temp += shannon(np.bincount(biopsy_Mutlist_temp)[x],cell_count_inBx[bx])
		SIBx_temp = float("{0:.3f}".format(SIBx_temp))
		SIBx.append(-SIBx_temp)
	for x in range (0, np.amax(aggregate_biopsy)):
		SIBx_agg += shannon(np.bincount(aggregate_biopsy)[x],np.sum(cell_count_inBx))
	SIBx_agg = float("{0:.3f}".format(SIBx_agg))
	return SIBx, -SIBx_agg

##################################################################

size = 100 #size of the array
time = 500
#pick X random points, then find all the other elements within a range, r, of the point
biopsy_num = 20 #desired number of biopsies
r = 10 #euclidean distance from random point that you include in biopsy
SI1 = 0 #placeholders for Shannon Index values
total_mut1 = np.zeros(size**2) #placeholders for mutation arrays

read_path = '../andrea_test/non-stem/text/'
write_path = '../andrea_test/non-stem/figs/'

################ GATHER AND PARSE DATA #################
SIBx = [] #list of shannon indices for each biopsy

#mutation flag data
data = open(read_path+'carriedMutation'+str(time)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('int')
CM1 = np.reshape(CA, (size,size))
N1 = np.count_nonzero(CA)
for species in range (1, np.amax(CA)): SI1 = SI1 + shannon(np.bincount(CA)[species],N1)
SItrunc = float("{0:.4f}".format(SI1))

################ INITIALIZE COUNTERS #################

#'biopsy' at random some circle of cells
biopsy_sites = [] #a list of the sites of biopsy - ordered pairs

### PLOT
rcParams['figure.figsize'] = 11,11

'''plot histogram unaggregated biopsies, radius #1'''

r1 = 5
biopsy_sites = gather_biopsies(biopsy_num,r1)
SIBx, SIBx_agg = do_biopsies_aggregate(size, biopsy_num, r1, CM1, biopsy_sites)
meanBx1 = np.mean(SIBx)
stdBx1 = np.std(SIBx)
plt.subplot(2,2,4)
# plt.hist(SIBx, SIBx_agg, 10, histtype = 'bar')
weights = np.ones_like(SIBx)/len(SIBx)
bins = np.linspace(0, np.max(SIBx), 100)
n, bins, patches = plt.hist(SIBx, 10, histtype = 'bar', weights = weights)
n, bins, patches = plt.hist(SIBx_agg, 10, histtype = 'bar', alpha = 0.7, color = 'r', \
	label = ('Aggregate, SI = '+str(SIBx_agg)))

plt.xlabel('Shannon Index')	
plt.ylabel('frequency')
plt.ylim([0, 1.2])
plt.legend(loc = 'upper center')
plt.title('S.I.s bxs of r = '+str(r1)+ ' \n mean:'+str(meanBx1)[:4]+' std:'+str(stdBx1)[:4])

plt.subplot(2,2,1)
ax1 = plt.pcolor(CM1, cmap='nipy_spectral', vmin = 0.001)
# plt.colorbar()
plt.title('Tumor size: '+str(N1)+' cells')
x,y = zip(*biopsy_sites)
# loop through all triplets of x-,y-coordinates and radius and plot a circle for each:
for x, y in zip(x, y):
    plt.gca().add_artist(Circle(xy=(x, y), radius=r1, alpha = 1, fill = False, color = 'w'))
plt.xlim([0, size])
plt.ylim([0, size])

'''plot histogram radius #2'''

r2 = 5
biopsy_sites = []
biopsy_sites = gather_biopsies(biopsy_num,r2)
SIBx = [] #list of shannon indices for each biopsy
SIBx, SIBx_agg = do_biopsies_aggregate(size, biopsy_num, r2, CM1, biopsy_sites)
meanBx2 = np.mean(SIBx)
stdBx2 = np.std(SIBx)

plt.subplot(2,2,3)
bins = np.linspace(0, np.max(SIBx), 100)
n, bins, patches = plt.hist(SIBx, 10, histtype = 'bar', weights = weights)
n, bins, patches = plt.hist(SIBx_agg, 10, histtype = 'bar', alpha = 0.7, color = 'r', \
	label = ('Aggregate, SI = '+str(SIBx_agg)))

plt.xlabel('Shannon Index')	
plt.ylabel('frequency')
plt.ylim([0, 1.2])
plt.legend(loc = 'upper center')
plt.title('S.I.s bxs of r = '+str(r2)+' \n mean:' +str(meanBx2)[:4]+' std:'+str(stdBx2)[:4])

'''PLOT associated CA and biopsy areas'''

plt.subplot(2,2,2)
ax1 = plt.pcolor(CM1, cmap='nipy_spectral', vmin = 0.001)
# plt.colorbar()
plt.title('Shannon Index: '+str(-SItrunc))
x,y = zip(*biopsy_sites)
for x, y in zip(x, y):
    plt.gca().add_artist(Circle(xy=(x, y), radius=r2, alpha = 1, fill = False, color = 'w'))
plt.xlim([0, size])
plt.ylim([0, size])

plt.savefig(write_path+'Indiv_vs_Aggreg Bx_N'+str(r1)+'.png', dpi = 500)
# plt.show()