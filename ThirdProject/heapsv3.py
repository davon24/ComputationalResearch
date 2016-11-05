#Code to implement Fibbonacci heaps
from numpy import *
import matplotlib.pyplot as plt
import random
import time as time

random.seed(3)
max_rank = 3

# Creates a class that will be used to define all nodes
class node:
    def __init__(self,key,children=[]):
        self.key = key
        self.children = children
   
    def rank(self):
    	return len(self.children)


# Function that combines two heaps
def combine(heap1, heap2):
	key1 = heap1.key
	key2 = heap2.key
# The root with the larger key is made a child of the other
	if key1 <= key2:
		heap1.children = heap1.children + [heap2]
		return heap1
	else:
		heap2.children = heap2.children + [heap1]
		return heap2


def arrange(heap,ranks):
	r = heap.rank()
# If the rank is higher than the max add one to the max
	try:
		ranks[r]
	except IndexError:
		ranks += [None]

# If another heap shares this rank combine the two
# and go through the process again
	if ranks[r] != None:
		new_heap = combine(heap,ranks[r])
		ranks[r] = None
		return arrange(new_heap, ranks)

# If nothing else shares the rank then take the spot
	else:
		ranks[r] = heap
		return ranks




def sorter(list_of_distances, max_rank):
	sorted_list = []
	t1 = time.time()
# Turns the distances into nodes
	heaps = [node(list_of_distances[i]) for i in range(len(list_of_distances))]
	t2 = time.time()
	print 'time to make heap: ' + str(t2-t1)


	ranked_list = [None for i in range(max_rank+1)] # List that heaps are stored in where the index is its rank
	for heap in heaps:
		arrange(heap,ranked_list)

	while len(sorted_list) < sqrt(len(list_of_distances)): # Continues until the first sqrt(N) are sorted
		min_key = 10**10
		min_index = 0

		# Finds the minimum
		for m in range(len(ranked_list)):
			if ranked_list[m] != None:
				keyi = ranked_list[m].key
				if keyi < min_key:
					min_key = keyi
					min_index = m


		childs = ranked_list[min_index].children # List of the minimum node's children
		sorted_list = sorted_list + [min_key] # Appends the minimum key to the sorted list
		ranked_list[min_index] = None	# Deletes minimum

		# The minimum's children become thier own heaps and are arranged with the others
		for kid in childs:
			arrange(kid, ranked_list)

		

	t3 = time.time()
	print "time to sort heap: " + str(t3-t2)
	return sorted_list


distanceList = [random.random() for i in range(100)]


t1 = time.time()
sorter(distanceList, 5)
t2 = time.time()


print "time to complete sort function: " + str(t2-t1)
