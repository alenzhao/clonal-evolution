''' Build phylogenies from parent,child lists

Jacob Scott 26 Dec 2015 '''


import numpy as np
from ete2 import Tree

data = open('output2c2.txt').read().replace(',',' ').replace('\n',' ')
x = data.split()
ParentChild = np.array(x).astype('str')
y = len(ParentChild)/3
ParentChild1 = np.reshape(ParentChild, (y,3))
# print ParentChild1

firsttwo = ParentChild1[:,0:2]
# print firsttwo

for row in range(0, len(firsttwo)): 
	for column in range(0,2): 
		firsttwo[row,column] = 'r'+firsttwo[row,column]
print firsttwo
# print type (firsttwo)

t = Tree() # Creates an empty tree
'''B = A.add_child(name="B")'''
r1 = t.add_child(name="r1")
lookup = {"r1": r1}

def sort_pairs(pair):
    # Extract integer after "r".
    return int(pair[0][1:])

for pair in sorted(firsttwo, key=sort_pairs):
    parentname = pair[0]
    childname = pair[1]
    if childname not in lookup:
        if parentname in lookup:
            # Add child.
            newchild = lookup[parentname].add_child(name = childname)
            lookup.update({childname: newchild})
        else:
            raise RuntimeError('Must not happen.')

print t.get_ascii(show_internal=True)
	


