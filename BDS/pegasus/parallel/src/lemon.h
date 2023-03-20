#ifndef LEMON_LIBRARY
#define LEMON_LIBRARY

#include <lemon/list_graph.h>
#include <lemon/gomory_hu.h>
#include <lemon/adaptors.h>
#include <lemon/connectivity.h>
#include <lemon/nauty_reader.h>

using namespace lemon;

typedef ListGraph::Node Node;
typedef ListGraph::Edge Edge;
typedef ListGraph::NodeIt NodeIt;
typedef ListGraph::EdgeIt EdgeIt;
// typedef ListGraph::NodeMap<int> NodeMap<int>;
// typedef ListGraph::NodeMap<double> NodeMap<double>;
// typedef ListGraph::EdgeMap<int> EdgeMap<int>;
// typedef ListGraph::EdgeMap<double> EdgeMap<double>;

ListGraph G; // Declare global Graph
// ListGraph::EdgeMap<int> cost(G); // Cost of the edges



#endif