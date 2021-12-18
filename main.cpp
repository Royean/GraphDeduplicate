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
#include "Timer.h"
#include "FPGrowth.h"

using namespace std;

extern bool compare(const NODE& x, const NODE& y); 

extern vector<string> split_line(string line, string delimiter);

extern int Rand(int i);

int main(int argc, char* args[]) {
    // select the input mode, bipartite or tripartite.
    bool isBipartite = atoi(args[1]) == 1 ? true : false;
    string filename(args[2]);
    printf("%s ", filename.c_str());
    vector<NODE> res;
    clock_t duration = 0;
    clock_t start = 0, end = 0;
    GraphDeduplicator* deduplicator = NULL;
    
    if(isBipartite){
        Bi2Tri bi2tri(filename);
        /*
        reduce the graph to smaller size.
        */
        start = clock();
        bi2tri.reduce();
        end = clock();
        duration += (end - start) /  CLOCKS_PER_SEC;
              
        bool needDeduplicate = true;
        if(needDeduplicate){
            /*
            compress the graph with duplication allowed.
            */
            start = clock();
            vector<NODE> tri = bi2tri.duplicate_allowed_compress();
            end = clock();
            duration += (end - start) /  CLOCKS_PER_SEC;
            deduplicator = new GraphDeduplicator(tri);
        }
        else{
            /*
            compress the graph without duplication.
            */
            start = clock();
            bi2tri.duplicate_free_compress();
            end = clock();
            duration += (end - start) /  CLOCKS_PER_SEC;
            return 0;
        }
    }
    else{
        /*
        adapt the input form for the deduplicator.
        */
        deduplicator = new GraphDeduplicator(filename);
    }

    exit(1);

    /**
     * @brief deduplicate 
     * 
     */
    start = clock();
    deduplicator->deduplicate();
    end = clock();
    duration += (end - start) /  CLOCKS_PER_SEC;
    
    // stringstream stream(args[1]);
    // float prob;
    // stream >> prob;
    // Bi2Tri bi2tri(1000, 2000, prob);
    // Bi2Tri bi2tri(prob, 4000, 5000);

    // test retail.dat

    // float minSupport;
    // FPGrowth *fp;   
    // minSupport = (float)40/100;
    // fp = new FPGrowth(minSupport);
    // fp->initFromFile(filename);
    // fp->outputToFile("retail_40.dat");
    // delete fp;

    // compress the graph directly, i.e., without duplication
    
    // printf("duplication time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);
    
    // start = clock();
    // // convert to a trigraph.
    // vector<NODE> res = bi2tri.convert();
    // end = clock();
    // printf("duplication time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);
    
    // clock_t end = clock();
    // printf("Input time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);
 
    // build the conflict graph from the three-column table here. 
    // The graph is stored in *neighbor* in the form of a 2-D array.
    // clock_t end = clock();

    // the input array is sorted since here, cannot reference the vertex with its original id. 
    // printf("Build graph time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);
    // duplicator.count_expand_edges();
    // cout << duplicator.LM  * 2 << " " <<duplicator.input.size()  << " "<< duplicator.exp_size << endl;
    // vector<int> mis = {0}, mvc;
    // for(int i = 1;i < duplicator.input.size();i++){
    //     mvc.push_back(i);
    // }

    // random shuffle the mvc
    // srand(unsigned(time(0)));
    // random_shuffle(mvc.begin(), mvc.end(), Rand);
    
    // report the expanded graph size.
    

    // naive deduplication version of dedup1 in sigmod
    // start = clock();
    // duplicator.dedup1();
    // end = clock();
    // printf("DeDup1 Deduplication time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);

    // greedy deduplication version of dedup1 in sigmod 
    // start = clock();
    // duplicator.greedyDedup();
    // end = clock();
    // printf("Greedy DeDup1 Deduplication time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);

    /* deduplicate the condensed graph */
    // start = clock();
    // sort(duplicator.input.begin(), duplicator.input.end(), compare); 
    // end = clock();
    // clock_t sort_time = end - start;
    // printf("Sort time:%fs\n", (double)(end - start) / CLOCKS_PER_SEC);
         
    // start = clock();
    // duplicator.build_conflict_graph();
    // duplicator.deduplicateBySearch(mvc, mis);
    // end = clock();
    // printf("Search Deduplication time:%fs\n", (double)(end - start + sort_time) / CLOCKS_PER_SEC);
    
    // start = clock();
    // duplicator.deduplicateBySetCover();
    // end = clock();
    // printf("Set Cover Deduplication time:%fs\n", (double)(end - start + sort_time) / CLOCKS_PER_SEC);  

    // start = clock();
    // duplicator.deduplicateByWeightedSetCover();
    // end = clock();
    // printf("Greedy Weighted Set Cover Deduplication time:%fs\n", (double)(end - start + sort_time) / CLOCKS_PER_SEC);  
    
    // duplicator.WeightedDeduplicateWithNegaEdge();

    // duplicator.deduplicateWithNegaEdge();
    return 0;
}