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
import os

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
		newpoint = [random.randint(r,size-r),random.randint(r,size-r)] #not including over the edge
		biopsy_sites.append(newpoint)
	return biopsy_sites

size = 100 #size of the array
time = 1000
#pick X random points, then find all the other elements within a range, r, of the point
biopsy_num = 3 #desired number of biopsies
r = 5 #euclidean distance from random point that you include in biopsy
SI1 = 0 #placeholders for Shannon Index values
total_mut1 = np.zeros(size**2) #placeholders for mutation arrays
detection_threshold = 0.6 #threshold for detection of clone/allele

#bit string data
data = open('../andrea_test/non-stem/text/genomes'+str(time)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('string')
Genomes = np.reshape(CA, (size,size))
genomelength = len(Genomes[0][0])
for entry in range(0, size**2):	total_mut1[entry] = np.array(sum_digits(CA[entry])).astype('int')
mut_array1 = np.reshape(total_mut1, (size,size))

#mutation flag data
data = open('../andrea_test/non-stem/text/carriedMutation'+str(time)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('int')
CM1 = np.reshape(CA, (size,size))
N1 = np.count_nonzero(CA)
for species in range (1, np.amax(CA)): SI1 = SI1 + shannon(np.bincount(CA)[species],N1)
SItrunc = float("{0:.4f}".format(SI1))

#'biopsy' at random some circle of cells
biopsy_sites = [] #a list of the sites of biopsy - ordered pairs
biopsied_cells = [] #a list of lists of biopsied cells
total_muts = np.zeros((biopsy_num,genomelength))
genomes_inBx = []
muts_inBx = []
total_mut_at_site = np.zeros((biopsy_num,genomelength))
muts_of_type = np.zeros((biopsy_num,genomelength))
positive_alleles = np.zeros(genomelength).astype('int') #to be used to store alleles that appear in ANY biopsy

point1 = [random.randint(r,size-r),random.randint(r,size-r)] #pick a random position at least r from the edge
biopsy_sites.append(point1)

biopsy_sites = gather_biopsies(biopsy_num,r)

# print biopsy_sites
cell_count_inBx = np.zeros(biopsy_num)

for bx in range(0, biopsy_num):
	biopsy_Genlist_temp = [] # for the bitstrings
	biopsy_Mutlist_temp = [] # for the mutation flags
	genomes_inBx_temp = Counter()
	muts_inBx_temp = Counter()
	punch = biopsy_sites[bx]
	for row in range(0, size):
		for column in range(0, size):
			a = (row,column)
			if distance.euclidean(a,punch) <= r: 
				biopsy_Genlist_temp.append(Genomes[column][row])
				biopsy_Mutlist_temp.append(CM1[column][row])
	for genome in biopsy_Genlist_temp:
		genomes_inBx_temp[genome]+=1
	for mutation in biopsy_Mutlist_temp:
		muts_inBx_temp[mutation]+=1
	cell_count_inBx[bx] = len(biopsy_Genlist_temp)

	genomes_inBx.append(genomes_inBx_temp)
	muts_inBx.append(muts_inBx_temp)

	for site in range(0, genomelength):  #iterate through each position in the genome
		site_list = []
		for cell in range(0, len(biopsy_Genlist_temp)): #over every cell in the biopsy
			cell_of_interest = biopsy_Genlist_temp[cell] #find all the position X genes for every cell
			if cell_of_interest[site] == '1' : positive_alleles[site] = 1 #flag alleles which appear in ANY biopsy
			site_list.append(cell_of_interest[site]) #make a list of them
		total_mut_at_site[bx][site] = sum_digits(site_list) #add up the total mutations at the site of interest
		if total_mut_at_site[bx][site]/len(biopsy_Genlist_temp) > detection_threshold: #find the percent positive
			total_muts[bx][site] = 1 #if > threshold, count as clonal
	for cell in range(0, len(biopsy_Mutlist_temp)): #over every cell in the biopsy
		for i in range(1,genomelength):
			if biopsy_Mutlist_temp[cell] == i: 
				muts_of_type[bx][i-1]+=1

np.savetxt('total_muts.txt', total_muts, fmt='%.0f')

allele_ID = np.linspace(1,genomelength,genomelength) #all possible alleles
truncation_list = [] #list of alleles that don't appear

#create truncated version of total_muts which only has entries at positions where positive_alleles = 1
for i in range (0,len(positive_alleles)):
	if positive_alleles[i] == 0:
		truncation_list.append(i)
total_muts_trunc = np.delete(total_muts, truncation_list, 1)
muts_of_type_trunc = np.delete(muts_of_type, truncation_list, 1)
total_mut_at_site_trunc = np.delete(total_mut_at_site, truncation_list, 1)
alleles_trunc = np.delete(allele_ID, truncation_list)

'''plot histograms'''

rcParams['figure.figsize'] = 20,10
for i in range(0,biopsy_num):
	plt.subplot(2, biopsy_num, i+1) 
	plt.bar(alleles_trunc, muts_of_type_trunc[i]/cell_count_inBx[i], align='center', alpha=0.4)
	plt.xticks(alleles_trunc, rotation = 315)
	plt.xlabel('frequecy of clone')	
	plt.xlabel('Unique mutation flag')
	plt.ylim([0, 1])
	plt.title('Biopsy # '+str(i+1)+': Position '+str(biopsy_sites[i]))

''' plot allele frequncies'''

for i in range(0,biopsy_num):
	colors = []
	allele_freq = total_mut_at_site_trunc[i]/cell_count_inBx[i]
	plt.subplot(2, biopsy_num, i+biopsy_num+1)
	for position in range(0,len(total_muts_trunc[i])): 
		if total_muts_trunc[i][position]==1: 
			colors.append('r')
		else: colors.append('b') #assign red color to alleles above detection threshold
	plt.bar(alleles_trunc, allele_freq, align='center', alpha=1, color=colors)
	plt.xticks(alleles_trunc, rotation = 315)
	plt.ylabel('frequency')	
	plt.xlabel('allele')

'''PLOT CA and biopsy areas'''

rcParams['figure.figsize'] = 13, 10
plt.figure()

ax1 = plt.pcolor(CM1, cmap='nipy_spectral', vmin = 0.001)
plt.colorbar()
plt.title('Shannon Index: '+str(-SItrunc))
x,y = zip(*biopsy_sites)

# loop through all triplets of x-,y-coordinates and radius and
# plot a circle for each:
for x, y in zip(x, y):
    plt.gca().add_artist(Circle(xy=(x, y), radius=r, alpha = 1, fill = False, color = 'w'))
plt.xlim([0, size])
plt.ylim([0, size])

plt.show()