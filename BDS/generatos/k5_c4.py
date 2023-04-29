#Generates a blowup pf the K5 in C4

def cycle(id):
	pre = 4*id
	return [ [pre, pre + 1, 1], [pre + 1, pre + 2, 1], [pre + 2, pre + 3, 1], [pre + 3, pre, 1]]

edges = []

for i in range(5):
	edges = edges + cycle(i)

edges.append([0, 4, 0])
edges.append([1, 8, 0])
edges.append([2, 12, 0])
edges.append([3, 16, 0])

edges.append([7, 9, 0])
edges.append([11, 13, 0])
edges.append([15, 17, 0])
edges.append([19, 5, 0])

edges.append([6, 14, 0])
edges.append([10, 18, 0])

print(20, len(edges))
for x in edges:
	print(x[0], x[1], x[2])