#ledder gedget 1

import sys

k = int(sys.argv[1])

def generate(id):
	edges = []

	pre = 8 * id;

	edges.append([pre + 5, pre + 6, 0])
	edges.append([pre + 5, pre + 7, 1])

	edges.append([pre + 6, pre + 8, 1])

	edges.append([pre + 7, pre + 8, 1])
	edges.append([pre + 7, pre + 9, 0])

	edges.append([pre + 8, pre + 10, 0])  
	
	edges.append([pre + 9, pre + 10, 1])
	edges.append([pre + 9, pre + 11, 1])

	edges.append([pre + 10, pre + 12, 1])

	edges.append([pre + 11, pre + 12, 0])  

	return edges

edges = []

edges.append([0, 1, 1])
edges.append([0, 2, 1])

edges.append([1, 2, 1])
edges.append([1, 3, 0])

edges.append([2, 4, 0])

edges.append([3, 4, 1])
edges.append([3, 5, 1])

edges.append([4, 6, 1])

for i in range(k):
	if i < k - 1 and i > 0:
		pre = i * 8 
		edges.append([pre + 11, pre + 13, 1])
		edges.append([pre + 12, pre + 14, 1])

	edges = edges + generate(i)	

pre = 4 + k * 8 

edges.append([pre - 1, pre + 1, 1])

edges.append([pre, pre + 2, 1])

edges.append([pre + 1, pre + 2, 1])
edges.append([pre + 1, pre + 3, 0])

edges.append([pre + 2, pre + 4, 0])

edges.append([pre + 3, pre + 4, 1])
edges.append([pre + 3, pre + 5, 1])

edges.append([pre + 4, pre + 5, 1])

edges.append([0, pre + 5, 0])

print(pre + 6, len(edges))
for x in edges:
	for j in x:
		print(j, end = ' ')
	print()	