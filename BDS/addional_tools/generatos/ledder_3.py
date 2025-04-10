#ledder gedget 1

import sys

k = int(sys.argv[1])

def generate(id):
	edges = []

	pre = 6 * id;

	edges.append([pre + 1, pre + 2, 1])
	edges.append([pre + 1, pre + 3, 0])

	edges.append([pre + 2, pre + 4, 0])

	edges.append([pre + 3, pre + 4, 1])
	edges.append([pre + 3, pre + 5, 1])

	edges.append([pre + 4, pre + 6, 1])

	edges.append([pre + 5, pre + 6, 0])

	return edges

edges = []

edges.append([0, 1, 1])
edges.append([0, 2, 1])

for i in range(k):
	if i < k - 1:
		pre = i * 6 
		edges.append([pre + 5, pre + 7, 1])
		edges.append([pre + 6, pre + 8, 1])

	edges = edges + generate(i)	

pre = (k - 1) * 6 

edges.append([pre + 5, pre + 7, 1])

edges.append([pre + 6, pre + 8, 1])

edges.append([pre + 7, pre + 8, 1])
edges.append([pre + 7, pre + 9, 0])

edges.append([pre + 8, pre + 10, 0])

edges.append([pre + 9, pre + 10, 1])
edges.append([pre + 9, pre + 11, 1])

edges.append([pre + 10, pre + 11, 1])


edges.append([0, pre + 11, 0])

print(pre + 12, len(edges))
for x in sorted(edges):
	for j in x:
		print(j, end = ' ')
	print()	