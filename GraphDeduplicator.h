#pragma once

// #include "Graph.h"
#include <vector>
#include <set>
#include <map>
#include <string>
using namespace std;


// #define BRUTE

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

        
        // store the degree 
        vector<int> left_deg;
        vector<int> right_deg;
        vector<int> virtual_left_deg;
        vector<int> virtual_right_deg;
        // status: inactive
        vector<int> left_status;
        vector<int> virtual_status;
        vector<int> right_status;
        // left-> virtual neighborhood
        vector<vector<int>> left_adj;
        vector<vector<int>> right_adj;
        vector<vector<int>> virtual_left_adj;
        vector<vector<int>> virtual_right_adj;

        GraphDeduplicator(string filename);
        GraphDeduplicator(vector<NODE> input);

        void deduplicate();
        
        vector<vector<int>> get_neighbour();
        vector<int> get_weight();
        pair<VVON, VVON> get_result();
        
        set<int> findCommonNeighbour(const set<int>& a, const set<int>& b);

        void buildReductionGraph();
        void degreeOneReduction();
        void twoHopReduction();

        void build_conflict_graph();
        void construct_graphTopology();
        void assign_vertex_weight(int k);
        // set<int> MVC(const vector<int>& weight, const vector<vector<int>>& neighbour, set<int>& mis);     

        void greedyDedup();
        void dedup1();
        
        void describe();
        void deduplicateBySearch(const vector<int>& mvc, const vector<int>& mis);
        void remove_edge_and_restore_connections(int target, bool left, set<int>& removed_edges, set<int>& common, set<pair<int,int>>& redundant_edges);
        void build_revertedIDX_for_realNodes();

        ui* Vnode = nullptr;
        // OUT_NODE* edge = new OUT_NODE[];
        vector<OUT_NODE> edges;
        vector<ui> next;

        void deduplicateBySetCover();
        void insert_edge(int u, int v, bool left);
       
        void deduplicateByWeightedSetCover();
        void deduplicateWithNegaEdge();
        void WeightedDeduplicateWithNegaEdge();

        // vector<int> getNeigbor(int k, method algo);
        void bfsTest(int iter);
        void k_core(int k);

        // void duplicateTest(int iterations);
        void degreeTest(int iter);
        void pagerankTest(int iter);

        void report_result();
        void print_graph();
        void count_expand_edges();
        void compute_mis_benefit(vector<int>& mis);
};

