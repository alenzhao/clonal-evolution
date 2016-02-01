''' Build phylogenies from parent,child lists

Jacob Scott 26 Dec 2015 '''


import numpy as np
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace, TextFace
from collections import Counter
from math import log as ln
import random as random

# read_path1 = '../../../../Thesis/phylogenies/experiment/non-stem/'
read_path = '../andrea_test/16M_trial/text/'
filename = 'output16M'

def layout(node):
    # if node.is_leaf():
    #     # Add node name to leaf nodes
    #     N = AttrFace("name", fsize=14, fgcolor="black")
    #     faces.add_face_to_node(N, node, 0)
    if "weight" in node.features:
        # Creates a sphere face whose size is proportional to node's
        # feature "weight"
        C = CircleFace(radius=node.weight, color="RoyalBlue", style="sphere")
        # Let's make the sphere transparent
        C.opacity = 0.6
        # And place as a float face over the tree
        faces.add_face_to_node(C, node, 0, position="float")

def sort_pairs(pair):
    # Extract integer after "r".
    return int(pair[0][1:])

'''PLAN to build a counter for all leaf nodes to make faces of appropriate size'''
# #mutation flag data
# data = open(read_path+'carriedMutation'+str(time)).read().replace(',','\n').replace('\n','')
# x = data.split()
# CA = np.array(x).astype('int')
# CM1 = np.reshape(CA, (size,size))

# #count up all members of each unique group
# for mutation in range(1, np.max(CM1)):
#     x,y = np.where(CM1 == mutation)

#life history information
data = open(read_path+filename+'.txt').read().replace(',',' ').replace('\n',' ')
x = data.split()
ParentChild = np.array(x).astype(str)
y = len(ParentChild)/5
ParentChild1 = np.reshape(ParentChild, (y,5))
firsttwo = ParentChild1[:,0:2] #chops off first line which encodes parameters of simulation and third column which is not yet used
parents = []
children = []

for row in range(0, len(firsttwo)): 
	for column in range(0,2): 
		firsttwo[row,column] = 'r'+firsttwo[row,column]

t = Tree() # Creates an empty tree

r1 = t.add_child(name="r1")
lookup = {"r1": r1}
prune_list = ['r1']

for pair in sorted(firsttwo, key=sort_pairs):
    parentname = pair[0]
    childname = pair[1]
    if childname not in lookup:
        if parentname in lookup:
            newchild = lookup[parentname].add_child(name = childname)
            lookup.update({childname: newchild})
            if parentname not in parents:
                prune_list.append(lookup[parentname])
            parents.append(parentname) #make list of unique terminal nodes (no children of children)
            children.append(newchild)
        else:
            raise RuntimeError('Must not happen.')

'''make a list of all leaves with no children, count them and then prune the tree of all leaves '''
prune_count = Counter(children) #counter than contains the number of children that each terminal node has]
# print(parents)
# print(prune_list)
# print(prune_count)
# for key in prune_count.keys():
#     # node_weight = ln(prune_count[key])
#     # node_weight = ln(n_children)
#     # print(node_weight,n_children)
#     # node = lookup[key]
#     node.add_features(weight=random.randint(50))

for n in t.traverse():
    n.add_face(TextFace(n.name, fsize = 16), column=0, position="branch-bottom")
    # n.add_features(weight=random.randint(0,20))

t.prune(prune_list)
# Create an empty TreeStyle
ts = TreeStyle()
# Set our custom layout function
ts.layout_fn = layout
# Draw a tree
# ts.mode = "c"  # this makes it circular
# False need to add node names manually
ts.show_leaf_name = False
ts.scale = 120

# Show branch data
ts.show_branch_length = False
ts.show_branch_support = True


# print (t.get_ascii(show_internal=True))

t.show(tree_style=ts)
