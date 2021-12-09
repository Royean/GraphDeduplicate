#pragma once
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include "Utility.h"

using namespace std;

class Bi2Tri{
    public:

    int e = 0;

    vector<vector<int> > left, right;

    Bi2Tri(string filename);
    Bi2Tri(int n, int m, float prob);

    Bi2Tri(float alpha, int n, int m);

    void compress();
    
    vector<NODE> convert(int redundancy);

    void summary();

};