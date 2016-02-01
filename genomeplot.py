'''
Plotting mutational load and allele frequency of entire tumor

Jacob Scott 21 Decemeber 2015

usage: python genomeplot.py read_path write_path filename size time

'''
import matplotlib.pyplot as plt
import numpy as np
from pylab import rcParams
import clonal_evolution_functions as cef
import sys

# Parse input arguments
read_path = str(sys.argv[1]) 
write_path = str(sys.argv[2]) 
filename = str(sys.argv[3])
size = int(sys.argv[4])
time = int(sys.argv[5])

my_cmap = plt.cm.get_cmap('nipy_spectral')
my_cmap.set_under('w')

# size = 1000 #size of the array
# time = 30000
total_mut1 = np.zeros(size**2)

# read_path = '../../../../Thesis/phylogenies/experiment/non-stem/text/'
# read_path = '../sweep_two/dot4/text/'
# write_path = '../figs/'
# filename = 'output_dot4_size1k_mf1em2_0Death_t30k'



data = open(read_path+filename+'.txt').read().replace(',',' ').replace('\n',' ')
x = data.split()
ParentChild = np.array(x).astype(str)
y = len(ParentChild)/5
ParentChild1 = np.reshape(ParentChild, (y,5))
firsttwo = np.array(ParentChild1[:,0:2]).astype(int) #chops off third-fifth which is not used here
parents = []
children = []

family_dict = {} #make dictionary of children and parents
for row in firsttwo:
	family_dict.update({row[1]: row[0]})

##########

#timepoint 1 clonal evolution
data = open(read_path+'carriedMutation'+str(time)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('int')
CM1 = np.reshape(CA, (size,size))
# print(CM1)

mutation_number = cef.total_mutation_map(CM1, size, family_dict)

rcParams['figure.figsize'] = 10,10

plt.subplot(2,2,2)
plt.pcolor(CM1, cmap='nipy_spectral', vmin = 0.001)
plt.title('ts '+str(time)+' - unique clones')
plt.xlabel('Unique clones')
plt.colorbar()
# plt.savefig(write_path+'UniqueMuts'+str(time)+'TEST.png', dpi = 500)

plt.subplot(2,2,1)
# plt.figure()
plt.pcolor(mutation_number, cmap='nipy_spectral', vmin = 0.001)
plt.title('ts '+str(time)+' - total muts')
plt.xlabel('total mutations')
plt.colorbar()

plt.subplot(2,2,3)
# plt.figure()
weightsTM = np.ones_like(mutation_number)/len(mutation_number)
binsTM = np.linspace(0, np.max(mutation_number), 100)
n, binsTM, patches = plt.hist(mutation_number, 10, histtype = 'bar', weights = weightsTM, normed = True)
plt.xlabel('Total mutations')
plt.xlim([1, np.amax(mutation_number)])

plt.subplot(2,2,4)
# plt.figure()
weights = np.ones_like(CM1)/len(CM1)
bins = np.linspace(0, np.max(CM1), 100)
n, bins, patches = plt.hist(CM1, 10, histtype = 'bar', weights = weights, normed = True)
# plt.hist(CM1, normed=True)
plt.xlabel('Unique clones')
plt.xlim([1, np.amax(CM1)])
# plt.yscale('log')

# plt.show()
plt.savefig(write_path+filename+'.png', dpi = 500)

# plt.show()