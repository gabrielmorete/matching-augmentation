n, m = map(int, input().split())
input()

edges = []
for i in range(m):
	a, b = map(int, input().split())
	edges.append([a, b])

input()
input()
s = input()
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

	val = []			

	t1, t2 = input().split('|')
	
	val.append(t1.split()[-1])

	cnt = 0
	lp = [float(x) for x in t2.split()]

	t1, t2 = input().split('|')
	
	cnt = 0
	ip = [float(x) for x in t2.split()]

	t1, t2 = input().split('|')
	
	cnt = 0
	bds = [float(x) for x in t2.split()]


	print(n, m)
	for i in range(m):
		print("{:10s} {:^5s} | {:^10s} {:^10s} {:^10s}".format( str(edges[i][0]) + " " + str(edges[i][1]), str(cost[i]), str(lp[i]), str(ip[i]), str(bds[i])))

	input() # blank line
	s = input()

