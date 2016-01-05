'''
Plotting multi-dimensional histogram shannon index for many values of biopsy radius

Jacob Scott 3 Jan 2016

'''
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import random
import pylab as P
from scipy.spatial import distance
# from collections import Counter
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

def gather_biopsies(biopsy_num, r):
	while len(biopsy_sites) < biopsy_num:
		# newpoint = [random.randint(r,size-r),random.randint(r,size-r)] #not including over the edge
		newpoint = [random.randint(0,size),random.randint(0,size)]#overlap is ok, near edge is ok
		biopsy_sites.append(newpoint)
	return biopsy_sites

def do_biopsies(size, biopsy_num, r, CM1, biopsy_sites): 
	cell_count_inBx = np.zeros(biopsy_num)
	for bx in range(0, biopsy_num):
		biopsy_Mutlist_temp = [] # for the mutation flags
		muts_inBx_temp = []
		punch = biopsy_sites[bx]
		for row in range(0, size):
			for column in range(0, size):
				a = (row,column)
				if distance.euclidean(a,punch) <= r: 
					biopsy_Mutlist_temp.append(CM1[column][row])
		cell_count_inBx[bx] = len(biopsy_Mutlist_temp)
		SIBx_temp = 0
		for x in range (0, np.amax(biopsy_Mutlist_temp)):
			SIBx_temp += shannon(np.bincount(biopsy_Mutlist_temp)[x],cell_count_inBx[bx])
		SIBx_temp = float("{0:.3f}".format(SIBx_temp))
		SIBx.append(-SIBx_temp)
	return SIBx

##################################################################

size = 100 #size of the array
time = 100
#pick X random points, then find all the other elements within a range, r, of the point
biopsy_num = 20 #desired number of biopsies
SI1 = 0 #placeholders for Shannon Index values
total_mut1 = np.zeros(size**2) #placeholders for mutation arrays
SIBx = [] #list of shannon indices for each biopsy
r_sweep = 3

read_path = '../andrea_test/non-stem/text/'
write_path = '../figs/multi_bx/'

################ GATHER AND PARSE DATA #################

data = open(read_path+'carriedMutation'+str(time)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype('int')
CM1 = np.reshape(CA, (size,size))
N1 = np.count_nonzero(CA)
for species in range (1, np.amax(CA)): SI1 = SI1 + shannon(np.bincount(CA)[species],N1)
SItrunc = float("{0:.4f}".format(SI1))

''' INITIALIZE COUNTERS '''

#'biopsy' at random some circle of cells
biopsy_sites = [] #a list of the sites of biopsy - ordered pairs
biopsied_cells = [] #a list of lists of biopsied cells
hist_list = [] #a list of histogram data

### PLOT

rcParams['figure.figsize'] = 8,8
z = np.zeros((r_sweep,biopsy_num))
for i in range(1,r_sweep+1):
	biopsy_sites = []
	SIBx = []
	biopsy_sites = gather_biopsies(biopsy_num,2*i)
	SIBx = do_biopsies(size, biopsy_num, 2*i, CM1, biopsy_sites)
	hist_list.append(SIBx)
data=[]
bins = np.linspace(0, np.max(hist_list), 100)
# for i in (np.arange(r_sweep)):
# 	data.append(hist_list[i])
# # data.append(hist_list[i] for i in (np.arange(r_sweep)))
# print data
# print hist_list[0]
# n, bins = P.hist([hist_list[0,:], hist_list[:,1], hist_list[:,2] ], 10, histtype='bar')



'''  plots separate plots for each histogram
for i in (np.arange(r_sweep)):
	plt.figure()
	plt.hist(hist_list[i], bins, alpha = 0.5, label = str(i))

'''

#makes a nice, but unreadable multicolored bar chart
x = np.linspace(1,r_sweep,r_sweep)
y = np.linspace(1,biopsy_num,biopsy_num)
x,y = np.meshgrid(x,y)

colors = ['r','g','b','k']
ci = 0

for i in range(0,r_sweep): 
	for j in range(0,biopsy_num): 
		z[i][j] = hist_list[i][j]
Z = np.reshape(z, (biopsy_num,r_sweep))
# print Z[:,1]
xs = np.arange(biopsy_num)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for y in (np.arange(r_sweep)):
	z = Z[:,y]
	ci += 1
	cz = ci%len(colors)
	ax.bar(xs, z, zs=y, zdir='y', alpha=0.8, color = colors[cz])

ax.set_xlabel('biopsies')
# ax.xticks(x, np.arange(biopsy_num))
ax.set_ylabel('radius of biopsy')
ax.set_zlabel('Shannon Index')


# plt.subplot(1,2,1)
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(x,y,Z, rstride=1, cstride=1, cmap='nipy_spectral', linewidth=0, antialiased=False)
# fig.colorbar(surf, shrink=0.5, aspect=5)
# plt.surf(SIBx)
# plt.xlabel('Shannon Index')	
# plt.ylabel('frequency')
# plt.title('S.I.s of r = '+str(r1)+ 'bx \n mean:'+str(meanBx1)+' std:'+str(stdBx1))

# plt.savefig(write_path+'NHist'+str(biopsy_num)+'Bx_r=2'+'thru_'+str(2*r_sweep)+'.png', dpi = 500)
plt.show()
