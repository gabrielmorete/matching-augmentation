import os
import subprocess

edges = []
adj = [[0] * 30] * 30

def run(extra):
	s = str(30) + ' ' + str(len(edges) + len(extra)) + "\n"
	for x in edges:
		s = s + str(x[0]) + ' ' + str(x[1]) + ' ' + str(x[2]) + '\n'
	for x in extra:
		s = s + str(x[0]) + ' ' + str(x[1]) + ' ' + str(x[2]) + '\n'

	# print(s)	
	command = "echo \"" + s + " \"  | ./main -stdio" # command to be executed

	out = str(subprocess.check_output(command, shell=True)).split('\\n')
	print(len(out))

	print(out[0])

	frac = int((out[-3].split()[2]))
	inte = int((out[-2].split()[2]))
	bds = (out[-2].split()[2])

	print(frac, inte, bds)




def addrec(v, extra, matched):
	if v == 30:
		run(extra)

	addrec(v + 1, adj, extra)

	if matched[v] == 0:
		for j in range(v + 1, 30):
			if adj[v][j] == 0 and matched[j] == 0:
				matched[v] = 1
				matched[j] = 1
				extra.append([v, j, 1])

				addrec(v + 1, extra, matched)
				
				extra.pop()
				matched[v] = 0
				matched[j] = 0



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

for i in range(len(edges)):
	edges[i][0] -= 1
	edges[i][1] -= 1


for x in edges:
	adj[x[0]][x[1]] = 1
	adj[x[1]][x[0]] = 1


extra = []
# addrec(1, adj, extra)


run(extra)