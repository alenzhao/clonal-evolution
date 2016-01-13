'''
Plotting population dynamics of entire tumor

Jacob Scott 13 Jan 2016

'''
import matplotlib.pyplot as plt
import numpy as np
from pylab import rcParams

read_path = '../andrea_test/stem/'
write_path = '../figs/andrea_flat/'

data = open(read_path+'output_test_counting.txt').read().replace(',',' ').replace('\n',' ')
x = data.split()
cellmatrix = np.array(x).astype(str)
y = len(cellmatrix)/5
data_shaped = np.reshape(cellmatrix, (y,5))
cells_over_time = data_shaped[:,2:].astype(int) #grabs stem, non-stem, and associated timesteps
total_over_time = np.zeros(len(cells_over_time))
for i in range(0,len(cells_over_time)):
	total_over_time[i] = (cells_over_time[i,0] + cells_over_time[i,1])
	print(cells_over_time[i,0],cells_over_time[i,1],cells_over_time[i,2], total_over_time[i])
plt.plot(cells_over_time[:,2],cells_over_time[:,1], color = 'k', linewidth = 5, label = 'non-stem')
plt.plot(cells_over_time[:,2],cells_over_time[:,0], color = 'b', linewidth = 5, label = 'stem')
plt.plot(cells_over_time[:,2],total_over_time, color = 'r', linewidth = 5, label = 'total')
plt.ylabel('Cells')
plt.xlabel('Time')
plt.legend(loc = 'upper left')

# plt.savefig(write_path+'Allele_freq_plot_NS_time'+str(time)+'D0_MF5em4_size1c.png', dpi = 500)
plt.show()