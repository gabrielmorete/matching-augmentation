#sebo vigen
import sys

def new_block(id):
	pre = (24) * (id - 1)

	edges = [];

	edges.append([pre + 2, pre + 3])
	edges.append([pre + 2, pre + 8])

	edges.append([pre + 3, pre + 4])
	edges.append([pre + 3, pre + 7])
	edges.append([pre + 3, pre + 8])

	edges.append([pre + 4, pre + 5])
	edges.append([pre + 4, pre + 6])

	edges.append([pre + 5, pre + 11])
	edges.append([pre + 5, pre + 17])

	edges.append([pre + 6, pre + 7])

	edges.append([pre + 7, pre + 9])

	edges.append([pre + 8, pre + 9])
	edges.append([pre + 8, pre + 13])
	edges.append([pre + 8, pre + 14])
	edges.append([pre + 8, pre + 20])

	edges.append([pre + 9, pre + 10])
	edges.append([pre + 9, pre + 12])

	edges.append([pre + 10, pre + 11])
	edges.append([pre + 10, pre + 12])

	edges.append([pre + 11, pre + 17])
	edges.append([pre + 12, pre + 13])

	edges.append([pre + 13, pre + 16])
	edges.append([pre + 13, pre + 17])

	edges.append([pre + 14, pre + 15])
	edges.append([pre + 14, pre + 20])

	edges.append([pre + 15, pre + 16])
	edges.append([pre + 19, pre + 19])

	edges.append([pre + 16, pre + 17])
	edges.append([pre + 16, pre + 18])

	edges.append([pre + 17, pre + 22])
	edges.append([pre + 17, pre + 23])

	edges.append([pre + 18, pre + 19])

	edges.append([pre + 19, pre + 21])

	edges.append([pre + 20, pre + 21])

	edges.append([pre + 21, pre + 22])

	edges.append([pre + 22, pre + 23])

	return edges


def new_connection(id):
	pre = (24) * (id - 1)

	edges = [];

	edges.append([pre + 20, pre + 26])
	edges.append([pre + 23, pre + 29])


	edges.append([pre + 21, pre + 24])
	edges.append([pre + 22, pre + 24])
	edges.append([pre + 24, pre + 25])


	edges.append([pre + 15, pre + 25])
	edges.append([pre + 25, pre + 28])
	edges.append([pre + 25, pre + 34])

	return edges


k = int(sys.argv[1])
opt = int(sys.argv[2])

edges = []

edges.append([0, 1])
edges.append([0, 2])
edges.append([0, 5])
edges.append([0, 10])

edges.append([1, 10])
edges.append([1, 4])

for i in range(1, k + 1):
	edges += new_block(i)

for i in range(1, k):
	edges += new_connection(i)

pre = (24) * (k - 1)

edges.append([pre + 15, pre + 21])

m = len(edges)
n = 0
for x in edges:
	n = max(n, x[0], x[1])

n += 1	

print(n, m)

if opt == 1:
	for x in edges:
		print(x[0], x[1], 1)	
else:
	for x in edges:
		print(x[0], x[1])