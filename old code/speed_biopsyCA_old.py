'''
Plotting CA, mutation map and mutational age.

For multiple biopsies plot allele frequency and detected alleles

detection_threshold tunable, as is biopsy_num and biopsy radius (r)

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

''' Take biopsies and return a list of the mutations present and number of cells '''
def return_biopsied_mutations(size, biopsy_num, r, CM1, biopsy_sites):
	area = 81#4*r**2
	biopsy_Mutlist = np.zeros((biopsy_num,area)).astype('int')
	cell_count_inBx = np.zeros(biopsy_num)
	for (row, col), cell in np.ndenumerate(CM1):
		a = (row,col)
		for bx in range(0, biopsy_num):
			punch = biopsy_sites[bx]
			if distance.euclidean(a,punch) <= r:
				biopsy_Mutlist[bx][cell_count_inBx[bx]] = cell
				cell_count_inBx[bx] += 1
	# print(biopsy_Mutlist)
	
	# biopsy_Mutlist = biopsy_Mutlist[biopsy_Mutlist>0]
	# print(biopsy_Mutlist)
	return biopsy_Mutlist, cell_count_inBx

''' Function to derive genome and count mutations in provided list of cells ''' 
def derive_genome_biopsy(biopsy_list, family_dict):
	derived_genomes_inBx = np.zeros(len(biopsy_list)).astype(str)
	# mutation_number_inBx = np.zeros((size,size)).astype(int)
	for position, cell in np.ndenumerate(biopsy_list):
		if cell == 0: continue
		temp_parent = 2
		bitstring = list('1')
		bitstring += (np.max(CM1)-1)*'0'
		if cell == 1:
			# derived_genomes_inBx[position] = bitstring
			derived_genomes_inBx[position] = ''.join(bitstring)
			# mutation_number_inBx[position] = sum_digits(bitstring)
			continue 
		else:
			while temp_parent > 1:
				temp_parent = family_dict[cell]
				bitstring[cell-1] = '1'
				if temp_parent == 1: break
				cell = family_dict[cell]
			# derived_genomes_inBx[position] = bitstring
			derived_genomes_inBx[position] = ''.join(bitstring)
			# mutation_number_inBx[position] = sum_digits(bitstring)
	return derived_genomes_inBx#, mutation_number_inBx

size = 50 #size of the array
time = 75
biopsy_num = 3 #desired number of biopsies
r = 5 #euclidean distance from random point that you include in biopsy
SI1 = 0 #placeholders for Shannon Index values
total_mut1 = np.zeros(size**2) #placeholders for mutation arrays
detection_threshold = 0.6 #threshold for detection of clone/allele
area = 4*r**2


read_path = '../andrea_test/non-stem/text/'
write_path = '../andrea_test/non-stem/figs/'
filename = 'output_speed_bx'

data = open('../andrea_test/non-stem/'+filename+'.txt').read().replace(',',' ').replace('\n',' ')
x = data.split()
ParentChild = np.array(x).astype(str)
y = len(ParentChild)/5
ParentChild1 = np.reshape(ParentChild, (y,5))
firsttwo = np.array(ParentChild1[:,0:2]).astype(int) #chops off third-fifth which is not used here

'''  function to derive genome and mutation number array from CM plot and family-history ''' 

family_dict = {} #make dictionary of children and parents
for row in firsttwo:
	family_dict.update({row[1]: row[0]})
# print(family_dict)

#mutation flag data
data = open(read_path+'carriedMutation'+str(time)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('int')
CM1 = np.reshape(CA, (size,size))
N1 = np.count_nonzero(CA)
for species in range (1, np.amax(CA)): SI1 = SI1 + shannon(np.bincount(CA)[species],N1)
SItrunc = float("{0:.4f}".format(SI1))

genomelength = np.max(CM1) #highest number unique mutation found in tumor

# counters etc. 
biopsy_sites = [] #a list of the sites of biopsy - ordered pairs
biopsied_cells = [] #a list of lists of biopsied cells
biopsy_Mutlist = np.zeros((biopsy_num,area)).astype('int')
cell_count_inBx = np.zeros(biopsy_num)

point1 = [random.randint(r,size-r),random.randint(r,size-r)] #pick a random position at least r from the edge
biopsy_sites.append(point1)

biopsy_sites = gather_biopsies(biopsy_num,r) #do biopies from function
biopsy_Mutlist, cell_count_inBx = return_biopsied_mutations(size, biopsy_num, r, CM1, biopsy_sites)#get mutations and cell numbers

total_muts = np.zeros((biopsy_num,genomelength))
total_mut_at_site = np.zeros((biopsy_num,genomelength))
muts_of_type = np.zeros((biopsy_num,genomelength))
positive_alleles = np.zeros(genomelength).astype('int') #to be used to store alleles that appear in ANY biopsy

''' test for mutations above threshold and write down what we find for each biopsy'''
for bx in range(0,biopsy_num):
	derived_genomes_inBx = derive_genome_biopsy(biopsy_Mutlist[bx], family_dict)
	for site in range(0, genomelength):  #iterate through each position in the genome
		site_list = []
		for genome in derived_genomes_inBx: #over every cell in the biopsy
			# print(len(genome))
			if genome[site] == '1' : positive_alleles[site] = 1 #flag alleles which appear in ANY biopsy
			site_list.append(genome[site]) #make a list of them
		total_mut_at_site[bx][site] = sum_digits(site_list) #add up the total mutations at the site of interest	
		if total_mut_at_site[bx][site]/len(derived_genomes_inBx) > detection_threshold: #find the percent positive
			total_muts[bx][site] = 1 #if > threshold, count as clonal
	for cell in range(0, len(derived_genomes_inBx)): #over every cell in the biopsy
		for i in range(1,genomelength):
			if derived_genomes_inBx[i] == 1: 
				muts_of_type[bx][i-1]+=1

np.savetxt('total_muts.txt', total_muts, fmt='%.0f')

allele_ID = np.linspace(1,genomelength,genomelength) #all possible alleles
truncation_list = [] #list of alleles that don't appear

#create truncated version of total_muts which only has entries at positions where there exists a positive_alleles = 1
for i in range (0,len(positive_alleles)):
	if positive_alleles[i] == 0:
		truncation_list.append(i)
total_muts_trunc = np.delete(total_muts, truncation_list, 1)
muts_of_type_trunc = np.delete(muts_of_type, truncation_list, 1)
total_mut_at_site_trunc = np.delete(total_mut_at_site, truncation_list, 1)
alleles_trunc = np.delete(allele_ID, truncation_list)



# '''plot histograms'''
# rcParams['figure.figsize'] = 10,10
# plt.hist(CM1, normed=True)


plt.figure()
rcParams['figure.figsize'] = 7,7
for i in range(0,biopsy_num):
	# plt.figure()
	plt.subplot(2, biopsy_num, i+1) 
	plt.bar(alleles_trunc, muts_of_type_trunc[i]/cell_count_inBx[i], align='center', alpha=0.4)
	plt.xticks(alleles_trunc, rotation = 315)
	plt.xlabel('frequecy of clone')	
	plt.xlabel('Unique mutation flag')
	plt.ylim([0, 1])
	plt.title('Biopsy # '+str(i+1)+': Position '+str(biopsy_sites[i]))

# for i in range(0,biopsy_num):
# 	alleles_truncNZ = np.nonzero(alleles_trunc)
# 	muts_of_type_truncNZ = np.nonzero(muts_of_type_trunc[i])
# 	print alleles_truncNZ
# 	print muts_of_type_truncNZ
# 	plt.subplot(2, biopsy_num, i+1) 
# 	plt.bar(alleles_truncNZ, muts_of_type_truncNZ/cell_count_inBx[i], align='center', alpha=0.4)
# 	plt.xticks(alleles_truncNZ, rotation = 315)
# 	plt.xlabel('frequecy of clone')	
# 	plt.xlabel('Unique mutation flag')
# 	plt.ylim([0, 1])
# 	plt.title('Biopsy # '+str(i+1)+': Position '+str(biopsy_sites[i]))

''' plot allele frequncies'''

##TODO USE np.nonzero() to display (within each for loop) only the non-zero elements per biopsy

for i in range(0,biopsy_num):
	# plt.figure()
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

ax1 = plt.pcolor(CM1.T, cmap='nipy_spectral', vmin = 0.001)
plt.colorbar()
plt.title('Shannon Index: '+str(-SItrunc))
x,y = zip(*biopsy_sites)

# loop through all triplets of x-,y-coordinates and radius and
# plot a circle for each:
for x, y in zip(x, y):
    plt.gca().add_artist(Circle(xy=(x, y), radius=r, alpha = 1, fill = False, color = 'w'))
plt.xlim([0, size])
plt.ylim([0, size])

plt.savefig(write_path+str(biopsy_num)+'TEST r= '+str(r)+'.png', dpi = 500)

# plt.show()