#include <iostream>
#include <vector>
#include <array>
#include <cassert>
#include "/Library/gurobi1001/macos_universal2/include/gurobi_c++.h"
#include <lemon/list_graph.h>

using namespace std;
using namespace lemon;  

typedef long long int ll;


/*
    Conventions
        - On lemon, graphs are 0 indexed;
*/






ListGraph G; // Declare global Graph
ListGraph::EdgeMap<int> cost(G); // Cost of the edges


// I will also save the input information
int _n, _m;
vector< array<int, 3> > _edges;



/*
    This function reads the Graph. The input format will be
    n m       number of nodes, edges
    a b c     edge bertween a, b with cost c
*/
void ReadInput(){
    int n, m;
    cin>>n>>m;
    
    _n = n;
    _m = m;

    for (int i = 0; i < n; i++){
        ListGraph::Node v = G.addNode();
        if (G.id(v) != i)
            cout<<"Error : node don't match id"<<endl;
        assert(G.id(v) == i);
    }

    for (int i = 0; i < m; i++){
        int a, b, c;
        cin>>a>>b>>c;

        a--; // 0-indexed
        b--; 

        // _edges.push_back({a, b, c});

        ListGraph::Edge e = G.addEdge(G.nodeFromId(a), G.nodeFromId(b));
        cost[e] = c;
    }

    // for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
    //     cout<<G.id(G.u(e))<<' '<<G.id(G.v(e))<<' '<<cost[e]<<endl;
}




// class IntegerCut: public GRBCallback {
//     public:
//         GRBVar** vars;
//         int n;
//         IntegerCut(GRBVar** xvars, int xn) {
//             vars = xvars;
//             n = xn;
//         }
//     protected:
//     void callback() {
//         try {
//             if (where == GRB_CB_MIPSOL) {
//                 double *x[n];
//                 for (i = 0; i < n; i++)
//                     x[i] = getSolution(vars[i], n);
                
//                 vector<int> vis(n, 0);
//                 dfs(0, n, x, vis);
             
//                 bool ok = 1;
//                 for (int i = 0; i < n; i++)
//                     if (vis[i] == 0)
//                         ok = 0;
//                 if (!ok) {
//                     GRBLinExpr expr = 0;
//                     for (int i = 0; i < n; i++)
//                         for (int j = 0; j < n; j++)
//                             expr += vars[i][j];
//                     addLazy(expr >= 2);
//                 }


//         } catch (GRBException e) {
//             cout << "Error number: " << e.getErrorCode() << endl;
//             cout << e.getMessage() << endl;
//         } catch (...) {
//             cout << "Error during callback" << endl;
//         }
//     }
// };






// // This function receivs a graph and returns a optmal fractional solution
// vector<double> IntegerSolution(Graph G){
//     try {
        
//         int n = G.order();

//         // Create an environment
//         GRBEnv env = GRBEnv(true);
// //        env.set("LogFile", "mip1.log");
//         env.start();
// `
//         // Create an empty model
//         GRBModel model = GRBModel(env);
//         model.set(GRB_IntParam_LazyConstraints, 1);
        
//         GRBVar **x = NULL;
//         x = new GRBVar*[n];
//         for (int i = 0; i < n; i++)
//             x[i] = new GRBVar[n];
        
//         for (int i = 0; i < n; i++) {
//             sort(G.adj[v].begin(), G.adj[v].end());
//             int p = 0, sz = G.adj[v].size();
//             for (int j = 0; j < i; j++) {
//                 while (p < sz and G.adj[v][p].first < j) p++;
//                 if (p < sz and G.adj[v][p].first == j)
//                     vars[i][j] = model.addVar(0.0, 1.0, G.adj[v][p].second, GRB_BINARY, "x_"+to_string(i)+"_"+to_string(j));
//                 else
//                     vars[i][j] = model.addVar(0.0, 0.0, 0, GRB_BINARY, "x_"+to_string(i)+"_"+to_string(j));
                
//                 var[j][i] = var[i][j];
//             }
//         }

//         for (i = 0; i < n; i++) {
//               GRBLinExpr expr = 0;
//               for (j = 0; j < n; j++)
//                 expr += vars[i][j];
//               model.addConstr(expr >= 2, "cut2_"+itos(i));
//         }
        
        
//         // Optimize model
//         model.optimize();

//         cout << x.get(GRB_StringAttr_VarName) << " "
//          << x.get(GRB_DoubleAttr_X) << endl;
//         cout << y.get(GRB_StringAttr_VarName) << " "
//          << y.get(GRB_DoubleAttr_X) << endl;
//         cout << z.get(GRB_StringAttr_VarName) << " "
//          << z.get(GRB_DoubleAttr_X) << endl;

//         cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

//     } catch(GRBException e) {
//         cout << "Error code = " << e.getErrorCode() << endl;
//         cout << e.getMessage() << endl;
//     } catch(...) {
//         cout << "Exception during optimization" << endl;
//     }

//     return 0;
// }




signed main(){
    ReadInput();
}
