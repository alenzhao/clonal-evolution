''' algorithm to make neightborhood of given size'''

import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt


r = 3
biopsy_site = [6,6]


#creating indices to euclidean neighborhood (for biopsy)
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
	return neighboring_cells

neighboring_cells = neighborhood(biopsy_site,r)
print(neighboring_cells)
print(len(neighboring_cells))

plt.scatter(*zip(*neighboring_cells))
plt.show()