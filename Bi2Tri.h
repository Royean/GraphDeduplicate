#pragma once
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include "Utility.h"

using namespace std;

class Bi2Tri{
    public:

    int ln = 0, rn = 0, e = 0;
    int max_degree = 0;

    vector<vector<int> > left, right;
    int* L_Deg, *R_Deg;
    int* L_valid, *R_valid;

    Bi2Tri(string filename);

    Bi2Tri(int n, int m, float prob);

    Bi2Tri(float alpha, int n, int m);

    // void compress();
    vector<NODE> duplicate_allowed_compress();
    void duplicate_free_compress();
    
    void degree_one_reduce();
    void upperbound_reduce();
    void reduce();

    void subgraph_reduce();
    vector<NODE> convert();
    void summary();

};