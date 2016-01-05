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
		newpoint = [random.randint(0,size),random.randint(0,size)]
		biopsy_sites.append(newpoint)
	return biopsy_sites

def do_biopsies(size, biopsy_num, r, CM1, biopsy_sites): 
	cell_count_inBx = np.zeros(biopsy_num)
	for bx in range(0, biopsy_num):
		biopsy_Mutlist_temp = [] # for the mutation flags
		muts_inBx_temp = []
		punch = biopsy_sites[bx]
		for row in range(0, size):
			for column in range(0, size):
				a = (row,column)
				if distance.euclidean(a,punch) <= r: 
					biopsy_Mutlist_temp.append(CM1[column][row])
					
		cell_count_inBx[bx] = len(biopsy_Mutlist_temp)
		SIBx_temp = 0
		for x in range (0, np.amax(biopsy_Mutlist_temp)):
			SIBx_temp += shannon(np.bincount(biopsy_Mutlist_temp)[x],cell_count_inBx[bx])
		SIBx_temp = float("{0:.3f}".format(SIBx_temp))
		SIBx.append(-SIBx_temp)
	return SIBx

##################################################################

size = 1000 #size of the array
time = 1000
#pick X random points, then find all the other elements within a range, r, of the point
biopsy_num = 100 #desired number of biopsies
r = 10 #euclidean distance from random point that you include in biopsy
SI1 = 0 #placeholders for Shannon Index values
total_mut1 = np.zeros(size**2) #placeholders for mutation arrays

read_path = '../../../../Thesis/phylogenies/experiment/non-stem/text/'
write_path = '../figs/andrea_flat/'

################ GATHER AND PARSE DATA #################
# #bit string data
# data = open(read_path+'genomes'+str(time)).read().replace(',','\n').replace('\n','')
# x = data.split()
# CA = np.array(x).astype('string')
# Genomes = np.reshape(CA, (size,size))
# genomelength = len(Genomes[0][0])
# for entry in range(0, size**2):	total_mut1[entry] = np.array(sum_digits(CA[entry])).astype('int')
# mut_array1 = np.reshape(total_mut1, (size,size))
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
biopsied_cells = [] #a list of lists of biopsied cells
# total_muts = np.zeros((biopsy_num,genomelength))
# genomes_inBx = []
muts_inBx = []
# total_mut_at_site = np.zeros((biopsy_num,genomelength))
# muts_of_type = np.zeros((biopsy_num,genomelength))

# ### make first point
# point1 = [random.randint(r,size-r),random.randint(r,size-r)] #pick a random position at least r from the edge
# biopsy_sites.append(point1)


### PLOT
rcParams['figure.figsize'] = 11,11

'''plot histogram radius #1'''

r1 = 30
biopsy_sites = gather_biopsies(biopsy_num,r1)
SIBx = do_biopsies(size, biopsy_num, r1, CM1, biopsy_sites)
meanBx1 = np.mean(SIBx)
stdBx1 = np.std(SIBx)
# skewBx1 = np.skew(SIBx)

plt.subplot(2,2,3)
plt.hist(SIBx)
plt.xlabel('Shannon Index')	
plt.ylabel('frequency')
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

r2 = 40
biopsy_sites = []
biopsy_sites = gather_biopsies(biopsy_num,r2)
SIBx = [] #list of shannon indices for each biopsy
SIBx = do_biopsies(size, biopsy_num, r2, CM1, biopsy_sites)
meanBx2 = np.mean(SIBx)
stdBx2 = np.std(SIBx)
# skewBx2 = np.skew(SIBx)

plt.subplot(2,2,4)
plt.hist(SIBx)
plt.xlabel('Shannon Index')	
plt.ylabel('frequency')
plt.title('S.I.s bxs of r = '+str(r2)+' \n mean:' +str(meanBx2)[:4]+' std:'+str(stdBx2)[:4])

'''PLOT associated CA and biopsy areas'''

plt.subplot(2,2,2)
ax1 = plt.pcolor(CM1, cmap='nipy_spectral', vmin = 0.001)
plt.colorbar()
plt.title('Shannon Index: '+str(-SItrunc))
x,y = zip(*biopsy_sites)
for x, y in zip(x, y):
    plt.gca().add_artist(Circle(xy=(x, y), radius=r2, alpha = 1, fill = False, color = 'w'))
plt.xlim([0, size])
plt.ylim([0, size])

plt.savefig(write_path+str(biopsy_num)+'Bx_N'+str(r1)+'_'+str(r2)+'.png', dpi = 500)
# plt.show()