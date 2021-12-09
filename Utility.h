#ifndef _UTILITY_H_
#define _UTILITY_H_

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <queue>
#include <set>

using namespace std;
#define _BITSET_ //use bit set to represent the adjacency matrix
#define _STATISTIC_
#define _KERNEL_
#define _RECOLOR_
//#define _BRACH_ON_COLOR_

//#define NDEBUG
#include <cassert>

#ifdef _BITSET_
#define set_bit(array, pos) (((array)[(pos)>>3]) |= (1<<((pos)&7)))
#define reverse_bit(array, pos) (((array)[(pos)>>3]) ^= (1<<((pos)&7)))
#define test_bit1(array, pos) (((array)[(pos)>>3])&(1<<((pos)&7)))
#define test_bit(array, pos) ((((array)[(pos)>>3])>>((pos)&7))&1)
#else
#define set_bit(array, pos) (((array)[pos]) = 1)
#define reverse_bit(array, pos) (((array)[pos]) = 1- ((array)[pos]))
#define test_bit(array, pos) ((array)[pos])
#endif

using ui = unsigned int; // vertex type
using ept = unsigned long; // edge pointer type; unsigned int can be used to process upto two billion undirected edges

typedef vector<int> VI;
typedef unsigned int ui;

struct NODE_COMMON_NEIGHBOUR {
    int node_id;
    set<int> left_common_neighbour;
    set<int> right_common_neighbour;
    NODE_COMMON_NEIGHBOUR(int id){
		this->node_id = id;
	}
    NODE_COMMON_NEIGHBOUR(int id, set<int> lcn, set<int> rcn){
		this->node_id = id;
		this->left_common_neighbour = lcn;
		this->right_common_neighbour = rcn;
	}
};

struct NODE {
    ui input_node_id;
    set<int> left;
    set<int> right;
	set<int> left_del;
	set<int> right_del;
    vector<NODE_COMMON_NEIGHBOUR> common_neighbour;
    NODE(int id) {
		this->input_node_id = id;
	}
    NODE(int id, set<int> left, set<int> right){
		this->input_node_id = id;
		this->left = left;
		this->right = right;
	}
};

struct EDGE {
    int u;
    int v;
    EDGE(int u, int v){
		this->u = u;
		this->v = v;
	}
};

struct OUT_NODE {
    int isVirtual; // 0 is false, 1 is true
    int node_id;
    OUT_NODE(int isVirtual, int node_id){
		this->isVirtual = isVirtual;
		this->node_id = node_id;
	}
};



typedef vector<vector<OUT_NODE>> VVON;

// int Rand(int i);

// enum method {search = 1, cover = 2};


// NODE_COMMON_NEIGHBOUR::NODE_COMMON_NEIGHBOUR(int id) {
//     this->node_id = id;
// }

// NODE_COMMON_NEIGHBOUR::NODE_COMMON_NEIGHBOUR(int id, set<int> lcn, set<int> rcn) {
//     this->node_id = id;
//     this->left_common_neighbour = lcn;
//     this->right_common_neighbour = rcn;
// }

// NODE::NODE(int id) {
//     this->input_node_id = id;
// }

// NODE::NODE(int id, set<int> left, set<int> right) {
//     this->input_node_id = id;
//     this->left = left;
//     this->right = right;
// }

// EDGE::EDGE(int u, int v) {
//     this->u = u;
//     this->v = v;
// }

// OUT_NODE::OUT_NODE(int isVirtual, int node_id) {
//     this->isVirtual = isVirtual;
//     this->node_id = node_id;
// }


#define EPSINON 1e-6
#define pb push_back
#define mp make_pair

extern bool compare(const NODE& x, const NODE& y); 
extern int Rand(int i);
extern vector<string> split_line(string line, string delimiter);


class Utility {
public:
	static FILE *open_file(const char *file_name, const char *mode) {
		FILE *f = fopen(file_name, mode);
		if(f == nullptr) {
			printf("Can not open file: %s\n", file_name);
			exit(1);
		}

		return f;
	}

	static std::string integer_to_string(long long number) {
		std::vector<ui> sequence;
		if(number == 0) sequence.pb(0);
		while(number > 0) {
			sequence.pb(number%1000);
			number /= 1000;
		}

		char buf[5];
		std::string res;
		for(ui i = sequence.size();i > 0;i --) {
			if(i == sequence.size()) sprintf(buf, "%u", sequence[i-1]);
			else sprintf(buf, ",%03u", sequence[i-1]);
			res += std::string(buf);
		}
		return res;
	}
};

#endif
