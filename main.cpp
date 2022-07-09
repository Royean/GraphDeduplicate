#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <algorithm>
#include "Graph.h"
#include "GraphDeduplicator.h"
#include "Bi2Tri.h"
#include "Utility.h"
#include "FPGrowth.h"

using namespace std;

extern bool compare(const NODE& x, const NODE& y); 

extern vector<string> split_line(string line, string delimiter);

extern int Rand(int i);

int main(int argc, char* args[]) {
    char* filename = args[1];

    // read the table 
    clock_t start = clock(); 
    GraphDeduplicator duplicator(filename);
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

    // greedy deduplication version of dedup1 in sigmod 
    start = clock();
    duplicator.greedyDedup();
    end = clock();
    printf("Greedy DeDup1 Deduplication time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);

    // deduplicate the condensed graph
    start = clock();
    sort(duplicator.input.begin(), duplicator.input.end(), compare); 
    end = clock();
    clock_t sort_time = end - start;
    printf("Sort time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);


    start = clock();
    duplicator.build_conflict_graph();
    duplicator.deduplicateBySearch(mvc, mis);
    end = clock();
    printf("Search Deduplication time:%fs\n", (double)(end - start + sort_time) / CLOCKS_PER_SEC);
    
    // start = clock();
    // duplicator.deduplicateBySetCover();
    // end = clock();
    // printf("Set Cover Deduplication time:%fs\n", (double)(end - start + sort_time) / CLOCKS_PER_SEC);  
    
    duplicator.buildReductionGraph();

    #ifndef BRUTE
    start = clock();
    duplicator.degreeOneReduction();
    duplicator.twoHopReduction();
    end = clock();
    printf("Reduction time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);  
    #endif

    start = clock();
    duplicator.deduplicateByWeightedSetCover();
    end = clock();
    printf("Greedy Weighted Set Cover Deduplication time:%fs\n", (double)(end - start + sort_time) / CLOCKS_PER_SEC);  
    

    start = clock();
    duplicator.WeightedDeduplicateWithNegaEdge();
    end = clock();
    printf("Greedy Weighted Set Cover Deduplication time:%fs\n", (double)(end - start + sort_time) / CLOCKS_PER_SEC);  


    // float minSupport;
    // FPGrowth *fp;   
    // minSupport = (float)40/100;
    // fp = new FPGrowth(minSupport);
    // fp->initFromFile(filename);
    // fp->outputToFile("retail_40.dat");
    // delete fp;
    // printf("duplication time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC); 
    
    return 0;
}

