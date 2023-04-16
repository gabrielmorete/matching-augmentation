import os
import subprocess

edges = []
adj = [[0] * 30 for x in range(30)]

cnt = 0

gbds = 0;
gip = 0;

def run(extra):
	global cnt, gbds, gip 
	s = str(30) + ' ' + str(len(edges) + len(extra)) + "\n"
	for x in edges:
		s = s + str(x[0]) + ' ' + str(x[1]) + ' ' + str(x[2]) + '\n'
	for x in extra:
		s = s + str(x[0]) + ' ' + str(x[1]) + ' ' + str(x[2]) + '\n'


	command = "echo \"" + s + " \"  | ./main -stdio" # command to be executed

	out = str(subprocess.check_output(command, shell=True)).split('\\n')

	frac = int(out[-4].split()[2])
	ip = int(out[-3].split()[2])
	bds = int(out[-2].split()[2])

	cnt += 1

	if gbds < bds/frac:
		gbds = bds/frac

	if gip < ip/frac:
		gip = ip/frac
		
	if ip/frac > 1.33 or bds/frac >= 1.4:
		print("oi")
		ss = str(frac) + " " + str(ip) + " " + str(bds) + " " + str(ip/frac) + " " + str(bds/frac)+ "\n " + s
		command = "echo \"" + ss + " \" >> petersen_out" # command to be executed
		# print(command)
		subprocess.check_output(command, shell=True)

	if (cnt % 1000 == 0):
		command = "echo \"" + str(cnt) + " " + str(gip) + " " + str(gbds) + "\" > log"
		subprocess.check_output(command, shell=True)

			

# only add edges from the outer arc to the inner arc
def addrec(v, extra, matched):
	if v >= 15:
		run(extra)
		return

	addrec(v + 1, extra, matched)

	# print(matched)
	if matched[v] == 0:
		for j in range(15, 30):
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


for i in range(30):
	for j in range(30):
		adj[i][j] = 0

for x in edges:
	adj[x[0]][x[1]] = 1
	adj[x[1]][x[0]] = 1


# print(adj)

extra = []
matched = [0] * 30
addrec(0, extra, matched)

print(cnt)
# run(extra)