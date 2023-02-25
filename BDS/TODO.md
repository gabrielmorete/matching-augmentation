	* do not initialize all n^2 variables. Try to use vector.
	* make an array of expr to deg2 constraints
	* change dfs to stop using edge id and use flags
	* change tree_adj to use a adaptor
	* change parent array to be Node, initialize parent[v] = v, update child of root before calling the dfs.