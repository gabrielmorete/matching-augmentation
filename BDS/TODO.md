### TO DO List
* ~~do not initialize all n^2 variables.~~
* ~~make an array of expr to deg2 constraints~~
* ~~change dfs to stop using edge id and use flags~~
* ~~change tree_adj to use a adaptor~~
* ~~change parent array to be Node, initialize parent[v] = v~~
* ~~add compiler optimizations to the makefile~~
* ~~merge nauty file and solver~~
* ~~split in libraries~~
* ~~use tarjan algorithm as a MIP separator~~
* ~~add a call to last processed geng line~~
* ~~supress solver output~~
* ~~add inconsistency checker and counter (3 tries -> report to log)~~
* ~~check memory comsumption~~
	* ~~memory leaks on gurobi get sol(valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt)~~
* command line
	* ~~add start at~~
* optimize BDS
* ~~optimize solver~~
* ~~optimize memory~~
* ~~change model to just update objective~~
* ~~use BDS as a initial solution~~
	