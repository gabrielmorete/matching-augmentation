import sys

op = int(sys.argv[1])

n, m = map(int, input().split())
print(n, m)

for i in range(m):
	a, b = map (int, input().split())
	if op == 0:
		a -= 1
		b -= 1
	else:
		a += 1
		b += 1
	print(a, b)
	
					