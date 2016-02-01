import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from pylab import rcParams
import os
import clonal_evolution_functions as cef
from scipy.spatial import distance
import scipy.stats as stats
import sys

read_path = '../figs/patchiness'
write_path = '../figs/patchiness/'

data = open(read_path+'/dot2/'+'means.txt').read().replace('\n',' ')
means_dot2 = data.split()
means_dot2 = np.array(means_dot2).astype(float)

data = open(read_path+'/dot2/'+'stds.txt').read().replace('\n',' ')
stds_dot2 = data.split()
stds_dot2 = np.array(stds_dot2).astype(float)

data = open(read_path+'/dot2/'+'skews.txt').read().replace('\n',' ')
skews_dot2 = data.split()
skews_dot2 = np.array(skews_dot2).astype(float)


data = open(read_path+'/dot4/'+'means.txt').read().replace('\n',' ')
means_dot4 = data.split()
means_dot4 = np.array(means_dot4).astype(float)

data = open(read_path+'/dot4/'+'stds.txt').read().replace('\n',' ')
stds_dot4 = data.split()
stds_dot4 = np.array(stds_dot4).astype(float)

data = open(read_path+'/dot4/'+'skews.txt').read().replace('\n',' ')
skews_dot4 = data.split()
skews_dot4 = np.array(skews_dot4).astype(float)


data = open(read_path+'/dot6/'+'means.txt').read().replace('\n',' ')
means_dot6 = data.split()
means_dot6 = np.array(means_dot6).astype(float)

data = open(read_path+'/dot6/'+'stds.txt').read().replace('\n',' ')
stds_dot6 = data.split()
stds_dot6 = np.array(stds_dot6).astype(float)

data = open(read_path+'/dot6/'+'skews.txt').read().replace('\n',' ')
skews_dot6 = data.split()
skews_dot6 = np.array(skews_dot6).astype(float)


data = open(read_path+'/dot8/'+'means.txt').read().replace('\n',' ')
means_dot8 = data.split()
means_dot8 = np.array(means_dot8).astype(float)

data = open(read_path+'/dot8/'+'stds.txt').read().replace('\n',' ')
stds_dot8 = data.split()
stds_dot8 = np.array(stds_dot8).astype(float)

data = open(read_path+'/dot8/'+'skews.txt').read().replace('\n',' ')
skews_dot8 = data.split()
skews_dot8 = np.array(skews_dot8).astype(float)


data = open(read_path+'/one/'+'means.txt').read().replace('\n',' ')
means_one = data.split()
means_one = np.array(means_one).astype(float)

data = open(read_path+'/one/'+'stds.txt').read().replace('\n',' ')
stds_one = data.split()
stds_one = np.array(stds_one).astype(float)

data = open(read_path+'/one/'+'skews.txt').read().replace('\n',' ')
skews_one = data.split()
skews_one = np.array(skews_one).astype(float)



rcParams['figure.figsize'] = 10,20
data_labels = ['dot2', 'dot4', 'dot6', 'dot8', 'one']

plt.subplot(3,1,1)
plt.boxplot([means_dot2, means_dot4, means_dot6, means_dot8, means_one])
plt.xticks([1, 2, 3, 4, 5], data_labels)
plt.title('Means')

plt.subplot(3,1,2)
plt.boxplot([stds_dot2, stds_dot4, stds_dot6, stds_dot8, stds_one])
plt.xticks([1, 2, 3, 4, 5], data_labels)
plt.title('StdDevs')

plt.subplot(3,1,3)
plt.boxplot([skews_dot2, skews_dot4, skews_dot6, skews_dot8, skews_one])
plt.xticks([1, 2, 3, 4, 5], data_labels)
plt.title('Skews')

plt.savefig(write_path+'panMomentDist.png', dpi = 500)
