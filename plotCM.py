'''
Plotting CA, mutation map.

Jacob Scott 5 Jan 2016

'''
import matplotlib.pyplot as plt
import numpy as np
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

size = 100 #size of the array
t1 = 250
t2 = 500
SI1 = 0 #placeholders for Shannon Index values
SI2 = 0
total_mut1 = np.zeros(size**2) #placeholders for mutation arrays
total_mut2 = np.zeros(size**2)
my_cmap = plt.cm.get_cmap('nipy_spectral')
my_cmap.set_under('w')

read_path = '../andrea_test/non-stem/text/'
write_path = '../figs/'

#timepoint 1 clonal evolution
data = open(read_path+'carriedMutation'+str(t1)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('int')
CM1 = np.reshape(CA, (size,size))
N1 = np.count_nonzero(CA)
for species in range (1, np.amax(CA)): SI1 = SI1 + shannon(np.bincount(CA)[species],N1)
SI1trunc = float("{0:.4f}".format(SI1))

#timepoint 2 clonal evolution
data = open(read_path+'carriedMutation'+str(t2)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('int')
CM2 = np.reshape(CA, (size,size))
N2 = np.count_nonzero(CA)
for species in range (1, np.amax(CA)): SI2 = SI2 + shannon(np.bincount(CA)[species],N2)
SI2trunc = float("{0:.4f}".format(SI2))

""" PLOT """

rcParams['figure.figsize'] = 15,5

plt.subplot(1, 2, 1)
plt.pcolor(CM1, cmap='nipy_spectral', vmin = 0.001)
plt.title('S.I. = '+str(-SI1trunc)+', cells = '+str(N1))
plt.colorbar()

plt.subplot(1, 2, 2)
plt.pcolor(CM2, cmap='nipy_spectral', vmin = 0.001)
plt.title('S.I. = '+str(-SI2trunc)+', cells = '+str(N2))
plt.colorbar()

# plt.savefig(write_path+'CMplot_run1_'+str(t2)+'.png', dpi = 500)

plt.show()
