''' Build phylogenies from parent,child lists

Jacob Scott 26 Dec 2015 '''


import numpy as np
from ete2 import Tree, ClusterTree
from collections import Counter

read_path1 = '../../../../Thesis/phylogenies/experiment/non-stem/'
read_path = '../andrea_test/non-stem/'

data = open(read_path+'output2k_size500.txt').read().replace(',',' ').replace('\n',' ')
x = data.split()
ParentChild = np.array(x).astype('str')
y = len(ParentChild)/3
ParentChild1 = np.reshape(ParentChild, (y,3))

firsttwo = ParentChild1[:,0:2]
parents = []
prune_list = []

for row in range(0, len(firsttwo)): 
	for column in range(0,2): 
		firsttwo[row,column] = 'r'+firsttwo[row,column]
# print firsttwo
# print type(firsttwo)

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
            # print parentname
            newchild = lookup[parentname].add_child(name = childname)
            lookup.update({childname: newchild})
            if parentname not in parents:
                prune_list.append(parentname)
            parents.append(parentname) #make list of unique terminal nodes (no children of children)
        else:
            raise RuntimeError('Must not happen.')

'''make a list of all leaves with no children, count them and then prune the tree of all leaves '''
prune_count = Counter(parents) #counter than contains the number of children that each terminal node has
# print prune_count
# print prune_list
# prune_list = prune_list[1:]
# t.prune(prune_list)
# t.ladderize()
# print t.get_ascii(show_internal=True)

t.show()


