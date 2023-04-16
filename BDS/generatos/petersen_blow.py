edges = []

#triangles
for i in range(1, 31, 3):
	edges.append([i, i + 1, 1])
	edges.append([i, i + 2, 1])
	edges.append([i + 1, i + 2, 1])

# rays
for i in range(2, 15, 3):
	edges.append([i, i + 15, 0])	

# outer
edges.append([3, 4, 0])	
edges.append([6, 7, 0])	
edges.append([9, 10, 0])	
edges.append([12, 13, 0])	
edges.append([15, 1, 0])

#inner

edges.append([16, 22, 0])	
edges.append([24, 28, 0])	
edges.append([30, 19, 0])	
edges.append([21, 25, 0])	
edges.append([27, 18, 0])

print(30, len(edges))
for x in sorted(edges):
	for j in x:
		print(j, end = ' ')
	print()	
