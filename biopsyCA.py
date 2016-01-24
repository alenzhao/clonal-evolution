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

# parameters
size = 500 #size of the array
time = 500
biopsy_num = 3 #desired number of biopsies
r = 30 #euclidean distance from random point that you include in biopsy
SI1 = 0 #placeholders for Shannon Index values
total_mut1 = np.zeros(size**2) #placeholders for mutation arrays
detection_threshold = 0.7 #threshold for detection of clone/allele
# area = 4*r**2


read_path = '../andrea_test/non-stem/text/'
write_path = '../andrea_test/non-stem/figs/'
filename = 'outputSPEEDTEST'

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

#mutation flag data
data = open(read_path+'carriedMutation'+str(time)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('int')
CM1 = np.reshape(CA, (size,size))
N1 = np.count_nonzero(CA)
for species in range (1, np.amax(CA)): SI1 = SI1 + cef.shannon(np.bincount(CA)[species],N1)
SItrunc = float("{0:.4f}".format(SI1))

genomelength = np.max(CM1) #highest number unique mutation found in tumor

''' This section collects biopsies randomly of radius r, and returns a list of all the mutations within each'''
biopsy_sites = cef.gather_biopsies(biopsy_num, r, size) #do biopies
biopsied_mutations = []
for i in range(0,biopsy_num):
	neighborhood_temp = cef.neighborhood(biopsy_sites[i],r)
	biopsied_mutations.append(cef.return_bx_muts_from_sites(CM1, neighborhood_temp))

#counters needed
observed_muts = np.zeros((biopsy_num,genomelength)) 
total_mut_at_site = np.zeros((biopsy_num,genomelength)) # to count up the number of mutations at each site
muts_of_type = np.zeros((biopsy_num,genomelength)) # to count up the mutations of each type
positive_alleles = np.zeros(genomelength) #to be used to store alleles that appear in ANY biopsy
derived_genomes_inBx = []

''' test for mutations above threshold and write down what we find for each biopsy '''
for bx in range(0,biopsy_num):
	derived_genomes_inBx.append(cef.derive_genome_biopsy(biopsied_mutations[bx], family_dict, CM1))
	# print(derived_genomes_inBx[0][0][10])
	for site in range(0, genomelength):  #iterate through each position in the genome
		site_list = []
		for genome in derived_genomes_inBx[bx]: #over every cell in the biopsy
			if genome[site] == '1': 
				positive_alleles[site] = 1 #flag alleles which appear in ANY biopsy
			site_list.append(genome[site]) #make a list of them
		# print(total_mut_at_site[bx][site])
		total_mut_at_site[bx][site] = cef.sum_digits(site_list) #add up the total mutations at the site of interest	
		if total_mut_at_site[bx][site]/len(derived_genomes_inBx[0]) > detection_threshold: #find the percent positive
			observed_muts[bx][site] = 1 #if > threshold, count as clonal
	for i in range(0,len(biopsied_mutations[bx])): #over every cell in the biopsy
		if biopsied_mutations[bx][i] > 0: muts_of_type[bx][biopsied_mutations[bx][i]-1] += 1

# np.savetxt('observed_muts.txt', observed_muts, fmt='%.0f')

allele_ID = np.linspace(1,genomelength,genomelength) #all possible alleles
truncation_list = [] #list of alleles that don't appear

#create truncated version of observed_muts which only has entries at positions where there exists a positive_alleles = 1
for i in range (0,genomelength):
	if positive_alleles[i] == 0:
		truncation_list.append(i)
total_muts_trunc = np.delete(observed_muts, truncation_list, 1)
muts_of_type_trunc = np.delete(muts_of_type, truncation_list, 1)
total_mut_at_site_trunc = np.delete(total_mut_at_site, truncation_list, 1)
alleles_trunc = np.delete(allele_ID, truncation_list)

allele_ticks = np.linspace(1,len(alleles_trunc),len(alleles_trunc))
'''PLOTS'''
rcParams['figure.figsize'] = 15,7
plt.figure()

for i in range(0,biopsy_num):
	# plt.figure()
	plt.subplot(2, biopsy_num, i+1) 
	plt.bar(allele_ticks, muts_of_type_trunc[i]/len(derived_genomes_inBx[0]), align='center', alpha=0.4)
	# plt.setp(fig1, xticks = allele_ticks, xticklabels = alleles_trunc)
	plt.xticks(allele_ticks, rotation = 315)
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

#plot allele frequncies

##TODO USE np.nonzero() to display (within each for loop) only the non-zero elements per biopsy
for bx in range(0,biopsy_num):
	# plt.figure()
	colors = []
	allele_freq = total_mut_at_site_trunc[bx]/len(derived_genomes_inBx[0])
	plt.subplot(2, biopsy_num, bx+biopsy_num+1)
	for i in range(0,len(total_muts_trunc[bx])):
		if total_muts_trunc[bx][i] == 1: 
			colors.append('r') #assign red color to alleles above detection threshold
		else: colors.append('b') 
	plt.bar(allele_ticks, allele_freq, align='center', alpha=1, color=colors)
	# plt.setp(fig1, xticks = allele_ticks, xticklabels = alleles_trunc)
	plt.xticks(allele_ticks, rotation = 315)
	# plt.xticklabels(alleles_trunc)
	plt.ylabel('frequency')	
	plt.xlabel('allele')

#PLOT CA and biopsy areas

rcParams['figure.figsize'] = 13, 10
plt.figure()

ax1 = plt.pcolor(CM1.T, cmap='nipy_spectral', vmin = 0.001)
plt.colorbar()
plt.title('Shannon Index: '+str(-SItrunc))
x,y = zip(*biopsy_sites)

# loop through all triplets of x-,y-coordinates and radius and
# plot a circle for each:
for x, y in zip(x, y):
    plt.gca().add_artist(Circle(xy=(y, x), radius=r, alpha = 1, fill = False, color = 'w'))
plt.xlim([0, size])
plt.ylim([0, size])

# plt.savefig(write_path+str(biopsy_num)+'TEST2 r= '+str(r)+'.png', dpi = 500)

plt.show()
