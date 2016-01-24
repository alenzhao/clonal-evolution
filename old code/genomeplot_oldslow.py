'''
Plotting mutational load and allele frequency of entire tumor

Jacob Scott 21 Decemeber 2015

'''
import matplotlib.pyplot as plt
import numpy as np
from pylab import rcParams

my_cmap = plt.cm.get_cmap('nipy_spectral')
my_cmap.set_under('w')

size = 20 #size of the array
time = 25
total_mut1 = np.zeros(size**2)

def sum_digits(digit):
    return sum(int(x) for x in digit if x.isdigit())

# read_path = '../../../../Thesis/phylogenies/experiment/non-stem/text/'
read_path = '../andrea_test/non-stem/text/'
write_path = '../figs/andrea_flat/'

#bit string data
data = open(read_path+'genomes'+str(time)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype(str)
Genomes = np.reshape(CA, (size,size))
genomelength = len(Genomes[0][0])
for entry in range(0, size**2):	total_mut1[entry] = np.array(sum_digits(CA[entry])).astype('int')
mut_array1 = np.reshape(total_mut1, (size,size))

#timepoint 1 clonal evolution
data = open(read_path+'carriedMutation'+str(time)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('int')
CM1 = np.reshape(CA, (size,size))
print(CM1)

rcParams['figure.figsize'] = 10,10

plt.subplot(2,2,1)
# plt.figure()
plt.pcolor(mut_array1, cmap='nipy_spectral', vmin = 0.001)
plt.title('ts '+str(time)+' - mut distance')
plt.xlabel('Unique mutations')
plt.colorbar()

plt.subplot(2,2,2)
# plt.figure()
plt.pcolor(CM1, cmap='nipy_spectral', vmin = 0.001)
plt.title('ts '+str(time)+' - unique muts')
plt.xlabel('Unique mutations')
plt.colorbar()

plt.subplot(2,2,3)
# plt.figure()
weightsTM = np.ones_like(total_mut1)/len(total_mut1)
binsTM = np.linspace(0, np.max(total_mut1), 100)
n, binsTM, patches = plt.hist(total_mut1, 10, histtype = 'bar', weights = weightsTM, normed = True)
# plt.hist(total_mut1, normed=True)
plt.xlabel('Total mutations')
plt.xlim([1, np.amax(total_mut1)])
# plt.yscale('log')

plt.subplot(2,2,4)
# plt.figure()
weights = np.ones_like(CM1)/len(CM1)
bins = np.linspace(0, np.max(CM1), 100)
n, bins, patches = plt.hist(CM1, 10, histtype = 'bar', weights = weights, normed = True)
# plt.hist(CM1, normed=True)
plt.xlabel('Total mutations')
plt.xlim([1, np.amax(CM1)])
# plt.yscale('log')

plt.savefig(write_path+'Allele_freq_plot_NS_time'+str(time)+'D0_MF5em4_size1c.png', dpi = 500)

# plt.show()