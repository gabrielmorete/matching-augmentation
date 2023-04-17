import sys

opt = int(sys.argv[1])

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
	lp = [float(x) for x in t2.split()]

	t1, t2 = input().split('|')
	val.append(t1.split()[-1])
	ip = [int(x) for x in t2.split()]


	t1, t2 = input().split('|')
	val.append(t1.split()[-1])
	bds = [int(x) for x in t2.split()]


	if opt == 0:
		print(n, m)
		for i in range(m):
			print("{:10s} {:^5s} | {:^10f} {:^10d} {:^10d}".format( str(edges[i][0]) + " " + str(edges[i][1]), str(cost[i]), lp[i], ip[i], bds[i]))
		print(val)
	else:
		print(n, m)
		for i in range(m):
			print("{:s} {:s} {:f}".format( str(edges[i][0]) + " " + str(edges[i][1]), str(cost[i]), lp[i]))

	s = input()
