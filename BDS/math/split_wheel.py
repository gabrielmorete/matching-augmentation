# Python script to generate split wheel instances

import sys

n = int(sys.argv[1])
opt = int(sys.argv[2])

edges = []

if opt == 1:
	for i in range(1, n + 1):
		edges.append([0, 2 * i - 1])

	for i in range(1, 2 * n):
		if i % 2:
			edges.append([i, i + 1])
		else:
			edges.append([i, (i + 2) % (2 * n + 1)])

	edges.append([2 * n, 2])			

elif opt == 2:
	for i in range(1, n + 1):
		edges.append([0, 2 * i - 1, 1])

	for i in range(1, 2 * n):
		if i % 2:
			edges.append([i, i + 1, 0])
		else:
			edges.append([i, (i + 2) % (2 * n + 1), 1])

	edges.append([2 * n, 2, 1])

elif opt == 3:
	assert n % 2 == 0

	for i in range(1, n + 1):
		edges.append([0, 2 * i - 1])

	for i in range(1, 2 * n):
		if i % 2:
			edges.append([i, i + 1])
		else:
			edges.append([i, (i + 2) % (2 * n + 1)])

	edges.append([2 * n, 2])			


	for i in range(1, 2 * n, 4):
		edges.append([i, i + 2])
else:
	print("?")				

print(2 * n + 1, len(edges))

for i in range(len(edges)):
	for y in edges[i]:
		print(y, end = ' ')
	print()	

