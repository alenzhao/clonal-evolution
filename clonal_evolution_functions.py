''' Functions to be used for analysis of clonal evolution CA 

Jacob G Scott Jan 2016

'''
from math import log as ln
import random
import numpy as np
from scipy.spatial import distance

''' calculate the shannon index '''
def shannon(n, N):
    if n == 0:
        return 0
    else:
        return (float(n)/N) * ln(float(n)/N)

''' sum up the digits in a string '''
def sum_digits(digit):
    return sum(int(x) for x in digit if x.isdigit())

'''find neighbors of prescribed size for biopsy'''
def neighborhood(biopsy_site,r):
	neighboring_cells = []
	x = biopsy_site[0]
	y = biopsy_site[1]
	for i in range(0,r+1):
		for j in range(0,r+1):
			if distance.euclidean((x+i,y+j),biopsy_site) <= r and (x+i,y+j) not in neighboring_cells:
				neighboring_cells.append((x+i, y+j))
			if distance.euclidean((x-i,y+j),biopsy_site) <= r and (x-i,y+j) not in neighboring_cells:
				neighboring_cells.append((x-i, y+j))
			if distance.euclidean((x-i,y-j),biopsy_site) <= r and (x-i,y-j) not in neighboring_cells:
				neighboring_cells.append((x-i, y-j))
			if distance.euclidean((x+i,y-j),biopsy_site) <= r and (x+i,y-j) not in neighboring_cells:
				neighboring_cells.append((x+i, y-j))
	return neighboring_cells#returns ordered pairs

''' Take biopsies '''
def do_biopsies(size, biopsy_num, r, CM1, biopsy_sites):
	area = 4*r**2
	biopsy_Mutlist = np.zeros((biopsy_num,area)).astype('int')
	cell_count_inBx = np.zeros(biopsy_num)
	for row in range(0, size):
		for column in range(0, size):
			a = (row,column)
			for bx in range(0, biopsy_num):
				punch = biopsy_sites[bx]
				if distance.euclidean(a,punch) <= r:
					biopsy_Mutlist[bx][cell_count_inBx[bx]] = CM1[column][row]
					cell_count_inBx[bx] += 1
	for bx in range(0, biopsy_num):
		SIBx_temp = 0
		biopsy_Mutlist_temp = (biopsy_Mutlist[bx])[0:cell_count_inBx[bx]]
		for x in range (0, np.amax(biopsy_Mutlist_temp)):
			SIBx_temp += shannon(np.bincount(biopsy_Mutlist_temp)[x],cell_count_inBx[bx])
		SIBx_temp = float("{0:.3f}".format(SIBx_temp))
		SIBx.append(-SIBx_temp)
	return SIBx

''' Function to derive genome and count mutations in provided list of cells ''' 
def derive_genome_biopsy(biopsy_list, family_dict, CM1):
	# derived_genomes_inBx = np.zeros(len(biopsy_list)).astype('S300')
	derived_genomes_inBx = list(map(str, np.zeros(len(biopsy_list))))
	# print(derived_genomes_inBx[0])
	# mutation_number_inBx = np.zeros((size,size)).astype(int)
	# for position, cell in np.ndenumerate(biopsy_list):
	for biopsy in range(0,len(biopsy_list)):
		if biopsy_list[biopsy] == 0:
			bitstring = (np.max(CM1))*'0'
			derived_genomes_inBx[biopsy] = ''.join(bitstring)
			continue
		# temp_parent = 2
		bitstring = list('1')
		bitstring += (np.max(CM1)-1)*'0'
		if biopsy_list[biopsy] == 1:
			# derived_genomes_inBx[position] = np.array(bitstring)
			# print(derived_genomes_inBx[biopsy])
			derived_genomes_inBx[biopsy] = ''.join(bitstring)
			# mutation_number_inBx[position] = cef.sum_digits(bitstring)
			continue 
		else:
			temp_parent = family_dict[biopsy_list[biopsy]]
			bitstring[biopsy_list[biopsy]-1] = '1'
			while temp_parent > 1:
				temp_parent = family_dict[position]
				bitstring[temp_parent-1] = '1'
				if temp_parent == 1: break			
				# derived_genomes_inBx[position] = np.array(bitstring)
				# print(derived_genomes_inBx[position])
				# mutation_number_inBx[position] = cef.sum_digits(bitstring)
			derived_genomes_inBx[biopsy] = ''.join(bitstring)
	return derived_genomes_inBx#, mutation_number_inBx

''' Function to count mutations in entire array ''' 
def count_mutations(CM1, family_dict):
	mutation_number = np.zeros((size,size)).astype(int)
	for (row, col), cell in np.ndenumerate(CM1):
		if cell == 0: continue
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
		if cell == 0: continue
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


''' Identify centers of areas to biopsy and return them, r currently redundant '''
def gather_biopsies(biopsy_num, r, size):
	biopsy_sites = []
	point1 = [random.randint(r,size-1-r),random.randint(r,size-1-r)] #pick a random position at least r from the edge
	biopsy_sites.append(point1)
	while len(biopsy_sites) < biopsy_num:
		newpoint = [random.randint(r,size-1-r),random.randint(r,size-1-r)] #not including over the edge
		# newpoint = [random.randint(0,size),random.randint(0,size)] #overlap OK, over edge ok
		biopsy_sites.append(newpoint)
	return biopsy_sites


''' Take biopsies and return a list of the mutations present and number of cells '''
def return_bx_muts_from_sites(CM1, neighboring_cells):
	biopsy_Mutlist = np.zeros(len(neighboring_cells)).astype('int')
	for i in range (0,len(neighboring_cells)):
		pair = neighboring_cells[i]	
		biopsy_Mutlist[i] = CM1[pair[1]][pair[0]]
	return biopsy_Mutlist

''' Take biopsies and return a list of the mutations present and number of cells '''
def return_biopsied_mutations(size, biopsy_num, r, CM1, biopsy_sites):
	area = 4*r**2
	biopsy_Mutlist = np.zeros((biopsy_num,area)).astype('int')
	cell_count_inBx = np.zeros(biopsy_num)
	for (row, col), cell in np.ndenumerate(CM1):
		a = (row,col)
		for bx in range(0, biopsy_num):
			punch = biopsy_sites[bx]
			if distance.euclidean(a,punch) <= r:
				biopsy_Mutlist[bx][cell_count_inBx[bx]] = cell
				cell_count_inBx[bx] += 1
	return biopsy_Mutlist, cell_count_inBx


''' Gather biopsies from CA, aggregate all cells into one list and calculate SI ''' 
def do_biopsies_aggregate(size, biopsy_num, r, CM1, biopsy_sites):
	area = 4*r**2
	biopsy_Mutlist = np.zeros((biopsy_num,area)).astype('int')
	aggregate_biopsy = np.array([]).astype('int')
	cell_count_inBx = np.zeros(biopsy_num)
	SIBx_agg = 0
	for row in range(0, size):
		for column in range(0, size):
			a = (row,column)
			for bx in range(0, biopsy_num):
				punch = biopsy_sites[bx]
				if distance.euclidean(a,punch) <= r:
					biopsy_Mutlist[bx][cell_count_inBx[bx]] = CM1[column][row]
					cell_count_inBx[bx] += 1 
					aggregate_biopsy = np.append(aggregate_biopsy,CM1[column][row])
	for bx in range(0, biopsy_num):
		SIBx_temp = 0
		biopsy_Mutlist_temp = (biopsy_Mutlist[bx])[0:cell_count_inBx[bx]]
		for x in range (0, np.amax(biopsy_Mutlist_temp)):
			SIBx_temp += shannon(np.bincount(biopsy_Mutlist_temp)[x],cell_count_inBx[bx])
		SIBx_temp = float("{0:.3f}".format(SIBx_temp))
		SIBx.append(-SIBx_temp)
	for x in range (0, np.amax(aggregate_biopsy)):
		SIBx_agg += shannon(np.bincount(aggregate_biopsy)[x],np.sum(cell_count_inBx))
	SIBx_agg = float("{0:.3f}".format(SIBx_agg))
	return SIBx, -SIBx_agg