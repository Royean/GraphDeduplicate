#pragma once
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include "Utility.h"

using namespace std;

struct Edge{
    int dst;
    // int is_compressed;
    int uid;
    int state;
    Edge(int dst, int uid){
        this->dst = dst;
        this->uid = uid;   
        this->state = 1;
        // this->valid = 1;
        // this->is_compressed = 0;
    }
};

struct State{
    int is_compressed;

    State(){
        this->is_compressed = 0;
    }
};

class Bi2Tri{
    public:
    /**
     * @brief 
     * m is the edge size of undirected graph size.
     * n is the vertex number.
     * ln is the size of vertices on left side.
     * rn is the size of vertices on the right side.
     */
    // int n, m;
    int ln = 0, rn = 0;
    int max_degree = 0;
    // vector<int> start;
    // vector<int> Deg;

    int e = 0;
    vector<vector<Edge> > left, right;
    vector<State> undirected_edge_state;
    int* L_Deg, *R_Deg;
    // vector<vector<int>> left_compressed, right_compressed;

    int* L_valid, *R_valid;

    Bi2Tri(string filename);

    ~Bi2Tri();
    
    Bi2Tri(int n, int m, float prob);

    Bi2Tri(float alpha, int n, int m);

    // void compress();
    vector<NODE> duplicate_allowed_compress();
    vector<NODE> simplified_duplicate_allowed_compress();
    void duplicate_free_compress();
    
    void degree_one_reduce();
    void upperbound_reduce();
    void reduce();

    void edge_reduce();
    void buildCPGraph();
    // void subgraph_reduce();
    vector<NODE> convert(int redundancy);
    void summary();
};