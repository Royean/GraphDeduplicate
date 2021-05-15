#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include "Graph.h"
#include "GraphDeduplicator.h"
using namespace std;

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
    
    // report the expanded graph size.
    duplicator.count_expand_edges();

    // invoke the deduplication algorithm for comparison.
    // start = clock();
    // duplicator.dedup1();
    // end = clock();
    // printf("DeDup1 Deduplication time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);

    // // // deduplicate the condensed graph
    // start = clock();
    // duplicator.build_conflict_graph();
    // duplicator.deduplicateBySearch(mvc, mis);
    // end = clock();
    // printf("Search Deduplication time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);
    
    start = clock();
    duplicator.deduplicateBySetCover();
    end = clock();
    printf("Set Cover Deduplication time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);  

    // start = clock();
    // duplicator.deduplicateByAdvancedSearch();
    // end = clock();
    // printf("Set Cover Deduplication time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);  
    
    return 0;
}