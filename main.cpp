#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <algorithm>
// #include "Graph.h"
#include "GraphDeduplicator.h"

using namespace std;

bool compare(const NODE& x, const NODE& y) {
    set<int>::iterator xit = x.left.begin(), yit = y.left.begin(), xend = x.left.end(), yend = y.left.end();
    int x_min = *xit;
    int y_min = *yit;
    if (x_min < y_min) {
        return true;
    } 
    else if (x_min == y_min) {
        return *(--xend) < *(--yend); //  
    }
    return false;
}

int main(int argc, char* args[]) {
    // parse the parameters.
    char* filename = args[1];

    // read the table 
    clock_t start = clock(); 
    GraphDeduplicator duplicator(filename);
    // clock_t end = clock();
    // printf("Input time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);
 
    // build the conflict graph from the three-column table here. 
    // The graph is stored in *neighbor* in the form of a 2-D array.
    clock_t end = clock();

    // // the input array is sorted since here, cannot reference the vertex with its original id. 
    printf("Build graph time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);

    vector<int> mis = {0}, mvc;
    for(int i = 1;i < duplicator.input.size();i++){
        mvc.push_back(i);
    }
    // random shuffle the mvc
    srand(unsigned(time(0)));
    random_shuffle(mvc.begin(), mvc.end(), Rand);
    
    // report the expanded graph size.
    duplicator.count_expand_edges();

    // naive deduplication version of dedup1 in sigmod
    // start = clock();
    // duplicator.dedup1();
    // end = clock();
    // printf("DeDup1 Deduplication time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);

    // // greedy deduplication version of dedup1 in sigmod 
    // start = clock();
    // duplicator.greedyDedup();
    // end = clock();
    // printf("Greedy DeDup1 Deduplication time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);

    // deduplicate the condensed graph
    start = clock();
    sort(duplicator.input.begin(), duplicator.input.end(), compare); 
    end = clock();
    clock_t sort_time = end - start;
    printf("Sort time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);
         
    // start = clock();
    // duplicator.build_conflict_graph();
    // duplicator.deduplicateBySearch(mvc, mis);
    // end = clock();
    // printf("Search Deduplication time:%fs\n", (double)(end - start + sort_time) / CLOCKS_PER_SEC);
    
    // start = clock();
    // duplicator.deduplicateBySetCover();
    // end = clock();
    // printf("Set Cover Deduplication time:%fs\n", (double)(end - start + sort_time) / CLOCKS_PER_SEC);  

    start = clock();
    duplicator.deduplicateByWeightedSetCover();
    end = clock();
    printf("Greedy Weighted Set Cover Deduplication time:%fs\n", (double)(end - start + sort_time) / CLOCKS_PER_SEC);  
    

    duplicator.deduplicateWithNegaEdge();
    
    return 0;
}