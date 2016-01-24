'''
Deriving bit strings from carriedMutation map and mutation list

'''

import matplotlib.pyplot as plt
import numpy as np
from pylab import rcParams


my_cmap = plt.cm.get_cmap('nipy_spectral')
my_cmap.set_under('w')

size = 50 #size of the array
time = 75
total_mut1 = np.zeros(size**2)

def sum_digits(digit):
    return sum(int(x) for x in digit if x.isdigit())

# read_path = '../../../../Thesis/phylogenies/experiment/non-stem/text/'
read_path = '../andrea_test/non-stem/text/'
write_path = '../figs/andrea_flat/'
filename = 'output_speed_bx'

#bit string data
data = open(read_path+'genomes'+str(time)).read().replace(',','\n').replace('\n','')
x = data.split()
CA = np.array(x).astype(str)
Genomes = np.reshape(CA, (size,size))
genomelength = len(Genomes[0][0])
for entry in range(0, size**2):	total_mut1[entry] = np.array(sum_digits(CA[entry])).astype('int')
mut_array1 = np.reshape(total_mut1, (size,size))

# print(Genomes)
# print(mut_array1)

#timepoint 1 clonal evolution
data = open(read_path+'carriedMutation'+str(time)).read().replace(',','\n').replace('\n','')
x = data.split()
CM = np.array(x).astype('int')
CM1 = np.reshape(CM, (size,size))

data = open('../andrea_test/non-stem/'+filename+'.txt').read().replace(',',' ').replace('\n',' ')
x = data.split()
ParentChild = np.array(x).astype(str)
y = len(ParentChild)/5
ParentChild1 = np.reshape(ParentChild, (y,5))
firsttwo = np.array(ParentChild1[:,0:2]).astype(int) #chops off third-fifth which is not used here
parents = []
children = []

'''  function to derive genome and mutation number array from CM plot and family-history ''' 

family_dict = {} #make dictionary of children and parents
for row in firsttwo:
	family_dict.update({row[1]: row[0]})


''' Function to derive genome of entire array ''' 
def derive_genome_array(CM1, family_dict):
	derived_genome = np.zeros((size,size)).astype(str)
	for (row, col), cell in np.ndenumerate(CM1):
		temp_parent = 2
		bitstring = list('1')
		bitstring += (np.max(CM1)-1)*'0'
		if cell == 1:
			derived_genome[row][col] = ''.join(bitstring)
			continue 
		else:
			while temp_parent > 1:
				temp_parent = family_dict[cell]
				bitstring[cell-1] = '1'
				if temp_parent == 1: break
				cell = family_dict[cell]
			derived_genome[row][col] = ''.join(bitstring)
	return derived_genome

''' Function to derive genome and count mutations of provided list of cells ''' 
def derive_genome_biopsy(biopsy_list, family_dict):
	derived_genomes_inBx = np.zeros(len(biopsy_list)).astype(str)
	for position, cell in np.ndenumerate(biopsy_list):
		temp_parent = 2
		bitstring = list('1')
		bitstring += (np.max(CM1)-1)*'0'
		if cell == 1:
			derived_genomes_inBx[position] = ''.join(bitstring)
			continue 
		else:
			while temp_parent > 1:
				temp_parent = family_dict[cell]
				bitstring[cell-1] = '1'
				if temp_parent == 1: break
				cell = family_dict[cell]
			derived_genomes_inBx[position] = ''.join(bitstring)
	return derived_genomes_inBx

''' Function to count mutations in entire array ''' 
def count_mutations(CM1, family_dict):
	mutation_number = np.zeros((size,size)).astype(int)
	for (row, col), cell in np.ndenumerate(CM1):
		temp_parent = 2
		bitstring = list('1')
		bitstring += (np.max(CM1)-1)*'0'
		if cell == 1:
			mutation_number[row][col] = sum_digits(bitstring)
			continue 
		else:
			while temp_parent > 1:
				temp_parent = family_dict[cell]
				bitstring[cell-1] = '1'
				if temp_parent == 1: break
				cell = family_dict[cell]
			mutation_number[row][col] = sum_digits(bitstring)
	return mutation_number

''' Function to count mutations and derive genome in entire array ''' 
def count_derive_mutations(CM1, family_dict):
	mutation_number = np.zeros((size,size)).astype(int)
	derived_genomes = np.zeros(len(biopsy_list)).astype(str)
	for (row, col), cell in np.ndenumerate(CM1):
		temp_parent = 2
		bitstring = list('1')
		bitstring += (np.max(CM1)-1)*'0'
		if cell == 1:
			derived_genomes[row][col] = ''.join(bitstring)
			mutation_number[row][col] = sum_digits(bitstring)
			continue 
		else:
			while temp_parent > 1:
				temp_parent = family_dict[cell]
				bitstring[cell-1] = '1'
				if temp_parent == 1: break
				cell = family_dict[cell]
			derived_genomes[row][col] = ''.join(bitstring)
			mutation_number[row][col] = sum_digits(bitstring)
	return mutation_number,derived_genomes


derived_genome = derive_genome_array(CM1, family_dict)
print(derived_genome)


plt.pcolor(CM1, cmap='nipy_spectral', vmin = 0.001)
plt.xlabel('Unique mutations')
plt.colorbar()
plt.show()

# for row in range(0, len(firsttwo)): 
# 	for column in range(0,2): 
# 		firsttwo[row,column] = 'r'+firsttwo[row,column]

# for pair in sorted(firsttwo, key=sort_pairs):
#     parentname = pair[0]
#     childname = pair[1]
#     if childname not in lookup:
#         if parentname in lookup:
#             # print parentname
#             newchild = lookup[parentname].add_child(name = childname)
#             lookup.update({childname: newchild})
#             if parentname not in parents:
#                 prune_list.append(lookup[parentname])
#             parents.append(parentname) #make list of unique terminal nodes (no children of children)
#             children.append(newchild)
#         else:
#             raise RuntimeError('Must not happen.')

