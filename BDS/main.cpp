#include <iostream>
#include <vector>
#include "/Library/gurobi1001/macos_universal2/include/gurobi_c++.h"

using namespace std;

class Graph{
    private:
        int n, m;
    public:
        vector<vector<pair<int, int>>> adj;
        vector<array<int, 3>> edges;
        
    Graph() : n(0), m(0) {} Graph (int _n) n(_n), m(0) {adj.resize(n)}
    int order() { return n; }
    int size() { return m; }
    void add_edge(int a, int b, int c) {
        adj[a].push_back({b, c});
        adj[b].push_back({a, c});
        edg.push_back({c, a, b});
        m++;
    }
};

Graph ReadInput(){
    int n, m;
    cin>>n>>m;
    
    Graph G = Graph(n);
    
    for (int i = 0; i < m; i++){
        int a, b, c;
        cin>>a>>b>>c;
        add_edge(a, b, c);
    }
}




class IntegerCut: public GRBCallback {
    public:
        GRBVar** vars;
        int n;
        IntegerCut(GRBVar** xvars, int xn) {
            vars = xvars;
            n = xn;
        }
    protected:
    void callback() {
        try {
            if (where == GRB_CB_MIPSOL) {
                double *x[n];
                for (i = 0; i < n; i++)
                    x[i] = getSolution(vars[i], n);
                
                vector<int> vis(n, 0);
                dfs(0, n, x, vis);
             
                bool ok = 1;
                for (int i = 0; i < n; i++)
                    if (vis[i] == 0)
                        ok = 0;
                if (!ok) {
                    GRBLinExpr expr = 0;
                    for (int i = 0; i < n; i++)
                        for (int j = 0; j < n; j++)
                            expr += vars[i][j];
                    addLazy(expr >= 2);
                }


        } catch (GRBException e) {
            cout << "Error number: " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch (...) {
            cout << "Error during callback" << endl;
        }
    }
};






// This function receivs a graph and returns a optmal fractional solution
vector<double> IntegerSolution(Graph G){
    try {
        
        int n = G.order();

        // Create an environment
        GRBEnv env = GRBEnv(true);
//        env.set("LogFile", "mip1.log");
        env.start();
`
        // Create an empty model
        GRBModel model = GRBModel(env);
        model.set(GRB_IntParam_LazyConstraints, 1);
        
        GRBVar **x = NULL;
        x = new GRBVar*[n];
        for (int i = 0; i < n; i++)
            x[i] = new GRBVar[n];
        
        for (int i = 0; i < n; i++) {
            sort(G.adj[v].begin(), G.adj[v].end());
            int p = 0, sz = G.adj[v].size();
            for (int j = 0; j < i; j++) {
                while (p < sz and G.adj[v][p].first < j) p++;
                if (p < sz and G.adj[v][p].first == j)
                    vars[i][j] = model.addVar(0.0, 1.0, G.adj[v][p].second, GRB_BINARY, "x_"+to_string(i)+"_"+to_string(j));
                else
                    vars[i][j] = model.addVar(0.0, 0.0, 0, GRB_BINARY, "x_"+to_string(i)+"_"+to_string(j));
                
                var[j][i] = var[i][j];
            }
        }

        for (i = 0; i < n; i++) {
              GRBLinExpr expr = 0;
              for (j = 0; j < n; j++)
                expr += vars[i][j];
              model.addConstr(expr >= 2, "cut2_"+itos(i));
        }
        
        
        // Optimize model
        model.optimize();

        cout << x.get(GRB_StringAttr_VarName) << " "
         << x.get(GRB_DoubleAttr_X) << endl;
        cout << y.get(GRB_StringAttr_VarName) << " "
         << y.get(GRB_DoubleAttr_X) << endl;
        cout << z.get(GRB_StringAttr_VarName) << " "
         << z.get(GRB_DoubleAttr_X) << endl;

        cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    return 0;
}




signed main(){
    
}
