'''
Plotting CA, mutation map and mutational age.

Jacob Scott 23 Decemeber 2015

'''
import matplotlib.pyplot as plt
import numpy as np
import random
from scipy.spatial import distance
from collections import Counter
from matplotlib.patches import Circle

#sum up the digits in a string
def sum_digits(digit):
    return sum(int(x) for x in digit if x.isdigit())

size = 20 #size of the array
time = 200
#pick a X random points, then find all the other elements within a range, r, of the point
biopsy_num = 3 #desired number of biopsies
r = 3 #euclidean distance from random point that you include in biopsy

total_mut1 = np.zeros(size**2) #placeholders for mutation arrays

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

#'biopsy' at random some circle of cells
biopsy_sites = [] #a list of the sites of biopsy - ordered pairs
biopsied_cells = [] #a list of lists of biopsied cells
total_muts = np.zeros(genomelength)
genomes_inBx = []
muts_inBx = []
# genomes_inBx_temp = Counter()
# muts_inBx_temp = Counter()


point1 = [random.randint(r,size-r),random.randint(r,size-r)] #pick a random position at least r from the edge
biopsy_sites.append(point1)

while len(biopsy_sites) < biopsy_num:
	newpoint = [random.randint(r,size-r),random.randint(r,size-r)]
	#print newpoint
	distances = []
	for element in range(0, len(biopsy_sites)):
		distances.append(distance.euclidean(newpoint,biopsy_sites[element]))
	if min(distances) > 2*r: biopsy_sites.append(newpoint)

print biopsy_sites

for site in range(0, biopsy_num):
	biopsy_Genlist_temp = [] # for the bitstrings
	biopsy_Mutlist_temp = [] # for the mutation flags
	genomes_inBx_temp = Counter()
	muts_inBx_temp = Counter()
	punch = biopsy_sites[site]
	for row in range(0, size):
		for column in range(0, size):
			a = (row,column)
			if distance.euclidean(a,punch) < r: 
				biopsy_Genlist_temp.append(Genomes[column][row])
				biopsy_Mutlist_temp.append(CM1[column][row])
	for genome in biopsy_Genlist_temp:
		genomes_inBx_temp[genome]+=1
	for genome in biopsy_Mutlist_temp:
		muts_inBx_temp[genome]+=1

	genomes_inBx.append(genomes_inBx_temp)
	muts_inBx.append(muts_inBx_temp)

	for site in range(0, genomelength):  #iterate through each position in the genome
		site_list = []
		for cell in range(0, len(biopsy_Genlist_temp)): #over every cell in the biopsy
			cell_of_interest = biopsy_Genlist_temp[cell] #find all the position X genes for every cell
			site_list.append(cell_of_interest[site]) #make a list of them
		total_mut_at_site = sum_digits(site_list) #add up the total mutations at the site of interest
		if total_mut_at_site/len(biopsy_Genlist_temp) > 0.5: #find the percent positive
			total_muts[site] = 1 #if > threshold, count as clonal

'''PLOT STUFF'''

# collect the histograms of interest
mut_keys = muts_inBx[2].keys()
y_pos = np.arange(len(mut_keys))
performance = [muts_inBx[2][k] for k in mut_keys]
plt.subplot(2, 1, 2)
plt.bar(y_pos, performance, align='center', alpha=0.4)
plt.xticks(y_pos, mut_keys)
plt.xlabel('cells with mutation')
plt.xlabel('Unique mutation flag')
plt.title('histogram of mutations per biopsy')


plt.subplot(2, 1, 1)
# fig = plt.figure()
# ax1 = plt.subplot2grid((3,2), (0,0), colspan=2, rowspan = 2)
ax1 = plt.pcolor(CM1, cmap='nipy_spectral', vmin = 0.001)
plt.colorbar()

# initialize axis, important: set the aspect ratio to equal
plt.title('Total mutations')

x,y = zip(*biopsy_sites)

# loop through all triplets of x-,y-coordinates and radius and
# plot a circle for each:
for x, y in zip(x, y):
    plt.gca().add_artist(Circle(xy=(x, y), radius=r, alpha = 0.1))

plt.xlim([0, size])
plt.ylim([0, size])

plt.show()
