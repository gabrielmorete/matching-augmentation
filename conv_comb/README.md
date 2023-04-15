### Instructions to run the code

On pegasus, type on the terminal

$ module load gurobi

to load gurobi libraries. To compile the code use

$ make main

##### Usage
To run the algorithm, type on the terminal

$	./conv_comb frac_points_file int_points_file

Where 'frac_points_file' and 'int_points_file' are the names of the files.

##### Optional Flags

	frac_points_file int_points_file -verbose -strong -coef dividend divisor

	-strong
		Check strong conjecture. Only checks if the point has x(\delta(v)) = 2 for every v in V.
	
	-verbose
		Display the convex combinations.
	
	-coef a b
		The code finds the best coefficient. By setting the coef it will warn the user of the points such that the lowest coeeficient is above this threshold. By default the coef is 4/3.