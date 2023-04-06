# Python script to generate split wheel instances

import sys

n = int(sys.argv[1])
opt = int(sys.argv[2])

print(2 * n + 1, 3 * n)

if opt == 1:
	for i in range(1, n + 1):
		print(0, 2 * i - 1)

	for i in range(1, 2 * n):
		if i % 2:
			print(i, i + 1)
		else:
			print(i, (i + 2) % (2 * n + 1))

	print(2 * n, 2)			

elif opt == 2:
	for i in range(1, n + 1):
		print(0, 2 * i - 1, 1)

	for i in range(1, 2 * n):
		if i % 2:
			print(i, i + 1, 0)
		else:
			print(i, (i + 2) % (2 * n + 1), 1)

	print(2 * n, 2, 1)

elif opt == 3:
	assert n % 2 == 0

	for i in range(1, 2 * n, 4)
		print(i, i + 2)

else:
	print("?")				
