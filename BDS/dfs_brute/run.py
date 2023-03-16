import os
import subprocess
import string 

# os.listdir(path) lists files and directories
# subprocess.call("make dfs_brute", shell=True)
path = "."
files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
# files is a list with the name of all FILES
	
for name in files:
	if name[0] != 'g':
		continue

	g = open(name);
	n, m = map(int, g.readline().split())
	g.readline();

	edges = []
	for i in range(m):
		a, b = map(int, g.readline().split())
		edges.append([a, b]);

	g.readline();
	g.readline();
	s = g.readline();
	while s[0] != 'N':
		cost = [1] * m
		t1, t2 = s.split(':')
		m_id = t1

		for x in t2.split(','):
			a, b = map(int, x.strip().split())

			for i in range(m):
				if edges[i][0] == a and edges[i][1] == b:
					cost[i] = 0
				if edges[i][1] == a and edges[i][0] == b:
					cost[i] = 0
		
		t1, t2 = g.readline().split('|')
		

		cnt = 0
		lp = [float(x) for x in t2.split()]

		g.readline(); # integer sol	
		g.readline(); # BDS sol
		g.readline(); # blank line
		s = g.readline();

		# generate input
		file_in = str(n) + " " + str(m) + "\n"
		for i in range(m):
			file_in = file_in + str(edges[i][0]) + " " + str(edges[i][1]) +  " " + str(cost[1]) +  " " + str(lp[i]) + "\n"
		file_in = file_in + name + " " + m_id + "\n"	
		


		command ="echo \"" + file_in + " \"  | ./dfs_brute" # command to be executed
		out = str(subprocess.check_output(command, shell=True))

		print(name, m_id)






	# finished reading the graph


	# command ="./brute < " + path + x # command to be executed
	# out = str(subprocess.check_output(command, shell=True))
	# if len(out) > 3:
	# 	print("File : " + path + x + " is a counter exemple!!")
	# 	print(out)

