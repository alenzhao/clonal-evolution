'''
Plotting mutational load and allele frequency of entire tumor

Jacob Scott 21 Decemeber 2015

'''
import matplotlib.pyplot as plt
import numpy as np
from pylab import rcParams
import scipy
import csv

my_cmap = plt.cm.get_cmap('nipy_spectral')
my_cmap.set_under('w')

size = 1000 #size of the array
time = 3600
total_mut1 = np.zeros(size**2)

def sum_digits(digit):
    return sum(int(x) for x in digit if x.isdigit())

# read_path = '../../../../Thesis/phylogenies/experiment/non-stem/text/'
read_path = '../text/'
write_path = '../figures/'
filename = 'outputTEST'

mut_array1 = np.zeros((size,size))

# #bit string data speed up?
with open(read_path+'genomes'+str(time).strip(),'r') as dest_f:
    data_iter = csv.reader(dest_f,delimiter = ',', quotechar = '"')
    data = [data for data in data_iter]
CM1 = np.asarray(data, dtype = str)
CM1 = scipy.delete(CM1,-1, 1)
for row in range(0, size):	
	for column in range(0, size):
		mut_array1[row][column] = np.array(sum_digits(CM1[row][column])).astype('int')

#bit string data
# data = open(read_path+'genomes'+str(time)).read().replace(',','\n').replace('\n','')
# x = data.split()
# CA = np.array(x).astype(str)
# Genomes = np.reshape(CA, (size,size))
# genomelength = len(Genomes[0][0])
# for entry in range(0, size**2):	total_mut1[entry] = np.array(sum_digits(CA[entry])).astype('int')
# mut_array1 = np.reshape(total_mut1, (size,size))

#timepoint 1 clonal evolution
data = open(read_path+'carriedMutation'+str(time)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('int')
CM1 = np.reshape(CA, (size,size))

rcParams['figure.figsize'] = 10,10

plt.subplot(2,2,1)
# plt.figure()
plt.pcolor(mut_array1, cmap='nipy_spectral', vmin = 0.001)
plt.title('ts '+str(time)+' - mut distance')
plt.xlabel('Total mutations')
plt.colorbar()

plt.subplot(2,2,2)
# plt.figure()
plt.pcolor(CM1, cmap='nipy_spectral', vmin = 0.001)
plt.title('ts '+str(time)+' - unique muts')
plt.xlabel('Unique mutations')
plt.colorbar()

plt.subplot(2,2,3)
mut_array1 = np.reshape(mut_array1, (size**2))
weightsTM = np.ones_like(mut_array1)/(size**2)
binsTM = np.linspace(0, np.max(mut_array1), np.max(mut_array1))
n, binsTM, patches = plt.hist(mut_array1, 8, histtype = 'bar', normed = True, weights = weightsTM)
plt.xlabel('Total mutations')
plt.xlim([1, np.amax(mut_array1)])
# plt.ylim([0,1.1])
# weightsTM = np.ones_like(total_mut1)/len(total_mut1)
# binsTM = np.linspace(0, np.max(total_mut1), 100)
# n, binsTM, patches = plt.hist(total_mut1, 10, histtype = 'bar', weights = weightsTM, normed = True)
# plt.xlabel('Total mutations')
# plt.xlim([1, np.amax(total_mut1)])

plt.subplot(2,2,4)
CM1 = np.reshape(CM1, (size**2))
# plt.figure()
weights = np.ones_like(CM1)/(size**2)
bins = np.linspace(0, np.max(CM1), np.max(CM1))
n, bins, patches = plt.hist(CM1, 25, histtype = 'bar', normed = True, weights= weights)
# plt.hist(CM1, normed=True)
plt.xlabel('Unique mutations')
plt.xlim([1, np.amax(CM1)])
# plt.ylim([0,1.1])
# plt.yscale('log')

# plt.savefig(write_path+'Allele_freq_plot_'+str(filename)+'.png', dpi = 500)

plt.show()