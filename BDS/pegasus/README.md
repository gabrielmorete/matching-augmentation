### Instructions to run the code

On pegasus, type on the terminal

$ module load gurobi

to load gurobi libraries. To compile the code use

$ make main

##### Usage
	-stdio
		Read input on human readable format, check stdio_reader.cpp for detais
	-verbose
		Display extra information
	-log_start
		Read log_progress file and warmstart based on it.
		Can't be used with -start
		Warning: If running in parallel, log_progress may be empty if the code is interrupted.
	-all_matchings
		generates all matchings for a given graph
	-support
		only solves IP and BDS if the support of the LP solution is equal to the graph		
	-start n
		Start at the nth graph of geng output
		Can't be used with -log_start
	-thread t
		Paralelize the code using t threads (only for graph6 input).		

If run without the -stdio flag, the algorithm assumes graph6 input.

##### Logs
	The code generates 3 log files
	- log
		Found a feasible example
	- log_progress
		Keeps track of the best instances found until now
	- exception
		Keeps track of errors (will only be created if an error happens)		
