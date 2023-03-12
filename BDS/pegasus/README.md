### Instructions to run the code

On pegasus, type on the terminal

$ module load gurobi

to load gurobi libraries. To compile the code use

$ make main

##### Usage
	-stdio
		Read input on humam readable format, check stdio_reader.cpp for detais
	-verbose
		Display extra information
	-start n
		Start at the nth graph of geng output		

If run without the -stdio flag, the algorithm assumes graph6 input. Also, it brutes all matchings of the instance.