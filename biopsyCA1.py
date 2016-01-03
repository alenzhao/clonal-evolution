'''
Plotting CA, mutation map and mutational age.

Jacob Scott 23 Decemeber 2015

'''
import matplotlib.pyplot as plt
import numpy as np
import random
from scipy.spatial import distance

#sum up the digits in a string
def sum_digits(digit):
    return sum(int(x) for x in digit if x.isdigit())

size = 50 #size of the array
time = 500

'''
# choose subarray, if you like
x1 = 58
x2 = 60
y1 = 58
y2 = 60
'''

total_mut1 = np.zeros(size**2) #placeholders for mutation arrays

#load data
data = open('genomes'+str(time)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('string')
Genomes = np.reshape(CA, (size,size))
genomelength = len(Genomes[0][0])
for entry in range(0, size**2):	total_mut1[entry] = np.array(sum_digits(CA[entry])).astype('int')
mut_array1 = np.reshape(total_mut1, (size,size))
#print Genomes

#'biopsy' at random
biopsy_num = 3 #desired number of biopsies
biopsy_sites = [] #a list of the sites of biopsy - ordered pairs
biopsied_cells = [] #a list of lists of biopsied cells
total_muts = np.zeros(genomelength)



'''
#choose subarray
sample1 = Genomes[x1:x2,y1:y2]
'''

#instead, let's pick a random point, then find all the other elements within a range, r, of the point
r = 10 #euclidean distance from random point that you include in biopsyX


point1 = [random.randint(r,size-r),random.randint(r,size-r)] #pick a random position at least r from the edge
biopsy_sites.append(point1)
#print point1

while len(biopsy_sites) < biopsy_num:
	newpoint = [random.randint(r,size-r),random.randint(r,size-r)]
	#print newpoint
	distances = []
	for element in range(0, len(biopsy_sites)):
		distances.append(distance.euclidean(newpoint,biopsy_sites[element]))
	#print distances
	if min(distances) >= 2*r: biopsy_sites.append(newpoint)

print biopsy_sites

for site in range(0, biopsy_num):
	biopsy_list_temp = []
	punch = biopsy_sites[site]
	for row in range(0, size):
		for column in range(0, size):
			a = (row,column)
			if distance.euclidean(a,punch) <= r: 
				biopsy_list_temp.append(Genomes[row][column])
		#print biopsy_list_temp

	#biopsied_cells.append(biopsy_list_temp)
	for site in range(0, genomelength):  #iterate through each position in the genome
		site_list = []
		for cell in range(0, len(biopsy_list_temp)): #over every cell in the biopsy
			cell_of_interest = biopsy_list_temp[cell] #find all the position X genes for every cell
			site_list.append(cell_of_interest[site]) #make a list of them
			#print site_list
		total_mut_at_site = sum_digits(site_list) #add up the total mutations at the site of interest
		#print total_mut_at_site
		if total_mut_at_site/len(biopsy_list_temp) > 0.9: #find the percent positive
			total_muts[site] = 1 #if > threshold, count as clonal
	print total_muts

		
''' next step is to incorporate code for drawing biopsy sites from below into code to plot the Genomes
 so you can see a circle where the biopsy is taken'''


	


#plot
plt.pcolor(mut_array1, cmap='nipy_spectral', vmin = 0.001)
plt.title('Total mutations')
#plt.colorbar()

#print biopsy_sites
x ,y = zip(*biopsy_sites)
area = np.pi * (10*r)**2

#plt.Circle((25,25), radius=5, color='k', fill = 'False')
figure(figsize = (50,10))
plt.scatter(x,y, s=area, facecolors='none', edgecolors='k', linewidths=5)
plt.xlim([0, 50])
plt.ylim([0, 50])
# figure(figsize=(20,10))
plt.show()
'''
for row in range(0, size):
	for column in range(0, size):
		biopsy_list_temp = []
		a = (row,column)
		for site in range(0, biopsy_num):
			if distance.euclidean(a,biopsy_sites[site]) <= r: biopsy_list_temp.append(Genomes[row][column])
'''




'''

#timepoint 1
plt.subplot(2, 3, 1)
plt.pcolor(Cells1, cmap='nipy_spectral')
plt.title('Timestep '+str(t1)+' - cells: '+str(N1))
plt.colorbar()

#plt.savefig("images/3plot.png", dpi = 500)

plt.show()
'''