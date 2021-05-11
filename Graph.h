// 2021.3.2

#ifndef _GRAPH_H
#define _GRAPH_H

#include <vector>
#include <set>

#define INF 2147483647

using namespace std;

struct NODE_COMMON_NEIGHBOUR {
    int node_id;
    set<int> left_common_neighbour;
    set<int> right_common_neighbour;
    NODE_COMMON_NEIGHBOUR(int id);
    NODE_COMMON_NEIGHBOUR(int id, set<int> lcn, set<int> rcn);
};

struct NODE {
    int input_node_id;
    set<int> left;
    set<int> right;
	set<int> left_del;
	set<int> right_del;
    vector<NODE_COMMON_NEIGHBOUR> common_neighbour;
    NODE(int id);
    NODE(int id, set<int> left, set<int> right);
};


struct EDGE {
    int u;
    int v;
    EDGE(int u, int v);
};

struct OUT_NODE {
    int isVirtual; // 0 is false, 1 is true
    int node_id;
    OUT_NODE(int isVirtual, int node_id);
};

typedef vector<vector<OUT_NODE>> VVON;


// used for MIS

struct Vertex {
	Vertex() {}
	bool del = 0;
	int weight = 0;
	int neiSum = 0;
	int neiSize = 0;
	vector<int> adjacent;
	vector<char> adjacent_state;
};

typedef vector<int> VI;
typedef unsigned int ui;

class Graph {
public:
	int total_weight = 0;
	Graph();
	Graph(int vertexNum);
	~Graph();
	int n, m;

	void readGraph(vector<vector<int>>& graph, vector<int>& weight);
	void outputGraph();
	void htThree();
	void htThreeAll(vector<int>& mvc, vector<int>& mis);
	void dtThreeAll();
	void wtThreeAll();

	void htPvertex();
	void htDegreeOne();
	void htDegreeTwo();
	void htLowDeg();
	void htNeighborhood();
	void htCommon();
	void bfs();
	

	void reconstruct_graph();

	void deleteVertex(int id, int par, VI& pone, VI& dominate, VI& degree_one, VI& degree_two);
	void deleteVertex(int id, int par, VI& pone, VI& degree_one, VI& degree_two);
	void deleteVertex(int v, VI& pone,VI& dominate, VI& degree_one, VI& degree_two);
	void deleteVertex(int v, VI& pone, VI& degree_one, VI& degree_two);

	void deleteVertex(int u, int* tri,int* edges, bool* adj, bool* dominate, VI& pone, VI& degree_one, VI& degree_two, VI& dominated);
	bool isNei(int v1, int v2);
	int commonWeight(VI& a, VI& b);
	int commonWeight(VI& a, VI& b, VI& common);
	void findTwoNeibor(int v, int& n1, int& n2);

	int oneEdge(vector<int>& a);
	int twoEdge(vector<int>& a);
	int threeEdge(vector<int>& a);
	int edgeHandler(vector<int>& list, int len);
	
	void edge_reduction(bool* adj);
	void bfsReduction(VI& pone, VI& degree_one, VI& degree_two, bool* adj);
	void common_neighbor_reduction(VI& pone, VI& degree_one, VI& degree_two, vector<pair<int, int>>& backTrack, vector<vector<int>>& circles, bool* adj);
	void neighborhood_reduction(VI& pone, VI& degree_one, VI& degree_two, vector<pair<int, int>>& backTrack, vector<vector<int>>& circles, bool* adj);
	void symmetric_folding(VI& pone, VI& degree_one, VI& degree_two, vector<pair<int, int>>& backTrack, vector<vector<int>>& circles, bool* adj, int* AdjShare);
	void initial_reduction(VI& pone, VI& degree_one, VI& degree_two, vector<pair<int, int>>& backTrack, vector<vector<int>>& circles, bool* adj, int* AdjShare);
	void compute_triangles(int* tri, int* edges, int* AdjShare, bool* adj, bool* dominate, VI& pone, VI& degree_one, VI& degree_two, VI& dominated, vector<pair<int, int>>& backTrack, vector<vector<int>>& circles);
	void clearPone(VI& pone, VI& degree_one, VI& degree_two, vector<pair<int,int>>& backTrack, vector<vector<int>>& circles,VI& dominated, int* edges,int* tri, bool* dominate, bool* adj);
	void clearPone(bool* adj,VI& pone, VI& degree_one, VI& degree_two, vector<pair<int, int>>& backTrack, vector<vector<int>>& circles, bool key);
	void clearsingle(bool* adj, VI& pone, VI& degree_one, VI& degree_two);
	int path_reduction(int* dp, int start, int end);
	void path_recover(int* V, int num);
};

#endif