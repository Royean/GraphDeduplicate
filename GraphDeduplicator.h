#ifndef _BUILD_GRAPH_H
#define _BUILD_GRAPH_H

#include "Graph.h"
#include <vector>
#include <set>
#include <map>
#include <string>

using namespace std;

class GraphDeduplicator {
    
    public:
        int N = 0;
        int LM = 0;
        int RM = 0;
        int raw_edge_tot = 0;
        int exp_size = 0;

        // map<int, int> m;
        vector<NODE> input;
        // vector<EDGE> extra_edge;
        vector<vector<int>> neighbour;
        vector<int> weight;
        VVON out_left;
        VVON out_right;

        GraphDeduplicator(char* filename);

        vector<vector<int>> get_neighbour();
        vector<int> get_weight();
        pair<VVON, VVON> get_result();
        
        set<int> findCommonNeighbour(const set<int>& a, const set<int>& b);
        void build_conflict_graph();
        void construct_graphTopology();
        void assign_vertex_weight(int k);
        set<int> MVC(const vector<int>& weight, const vector<vector<int>>& neighbour, set<int>& mis);

        void dedup1();
        
        void deduplicateBySearch(const vector<int>& mvc, const vector<int>& mis);
        void remove_edge_and_restore_connections(int target, bool left, set<int>& removed_edges, set<int>& common, set<pair<int,int>>& redundant_edges);
        void build_revertedIDX_for_realNodes();

        ui* Vnode = nullptr;
        // OUT_NODE* edge = new OUT_NODE[];
        vector<OUT_NODE> edges;
        vector<ui> next;
        void deduplicateBySetCover();
        void insert_edge(int u, int v, bool left);
       
        void report_result();
        void print_graph();
        int count_expand_edges();
        void compute_mis_benefit(vector<int>& mis);
};

#endif
