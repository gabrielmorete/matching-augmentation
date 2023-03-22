import os
import subprocess
import string 

# os.listdir(path) lists files and directories
# subprocess.call("make dfs_brute", shell=True)

path = "."
# files is a list with the name of all FILES

# variables for smallest example with gap >= 1.4

#variables for largest gap
name_best = "none"
match_best = "none"
gap_best = 1;
size_best = 10000

cnt = 0
for name in os.listdir(path):
	if os.path.isfile(os.path.join(path, name)) == 0 or name[0] != 'g':
		continue
		
	g = open(name)
	n, m = map(int, g.readline().split())
	g.readline()

	edges = []
	for i in range(m):
		a, b = map(int, g.readline().split())
		edges.append([a, b])

	g.readline()
	g.readline()
	s = g.readline()
	while len(s) > 0 and s[0] != 'N':
		cost = [1] * m
		t1, t2 = s.split(':')
		m_id = t1

		for x in t2.split(','):
			a, b = map(int, x.strip().split())
			for i in range(m):
				if (edges[i][0] == a) and (edges[i][1] == b):
					cost[i] = 0
				if (edges[i][1] == a) and (edges[i][0]) == b:
					cost[i] = 0

		t1, t2 = g.readline().split('|')
		
		cnt = 0
		lp = [float(x) for x in t2.split()]

		g.readline() # integer sol	
		g.readline() # BDS sol
		g.readline() # blank line
		s = g.readline()

		# generate input
		file_in = str(n) + " " + str(m) + "\n"
		for i in range(m):
			file_in = file_in + str(edges[i][0]) + " " + str(edges[i][1]) +  " " + str(cost[i]) +  " " + str(lp[i]) + "\n"
		file_in = file_in + name + " " + m_id + "\n"	
		
		# print("---")
		# print(file_in)

		command ="echo \"" + file_in + " \"  | ./dfs_brute" # command to be executed
		out = str(subprocess.check_output(command, shell=True))

		f_min, f_max, lp_val =  out.split()
		f_min = float(f_min[2:])
		f_max = float(f_max)
		lp_val = float(lp_val[:-3])

		# print(name, m_id, f_min, f_max, lp_val)

		if f_max > gap_best:
			gap_best = f_max
			name_best = name
			match_best = m_id
			size_best = m
		if f_max == gap_best and m < size_best:
			name_best = name
			match_best = m_id
			size_best = m
	cnt = cnt + 1		
	if cnt % 100000 == 0:
		print("Progress: ", cnt, " | ", gap_best,"|", name_best,"-",match_best)	

print("Largest gap:", gap_best,"|", name_best,"-",match_best)