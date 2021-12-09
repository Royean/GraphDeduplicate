#pragma once
#include <cstring>
#include <fstream>
#include <cmath>
#include <set>
#include <map>
#include <vector>
#include <queue>
#include <algorithm>
#include <iterator>
#include <utility>
#include <iostream>
#include <ctime>
#include "Utility.h"
#include "GraphDeduplicator.h"
#include "assert.h"


#define DEBUG
// using namespace std;

extern bool compare(const NODE& x, const NODE& y); 

extern vector<string> split_line(string line, string delimiter);

extern int Rand(int i);

inline set<int> set_intersec(set<int>& a, set<int>& b) {
    set<int> tmp;
    set_intersection(a.begin(), a.end(), b.begin(), b.end(), insert_iterator<set<int>> (tmp, tmp.begin()));
    return tmp;
}

inline set<int> set_minus(set<int>& parent, set<int>& removed) {
    set<int> temp;
    set_difference(parent.begin(), parent.end(), removed.begin(), removed.end(), insert_iterator<set<int>> (temp, temp.begin()));
    return temp;
}

void print_pair(pair<int, int> t) {
    cout << "(" << t.first << "," << t.second << ")" << endl;
}

template<class ForwardIterator>
inline size_t argmin(ForwardIterator first, ForwardIterator last)
{
    return std::distance(first, std::min_element(first, last));
}

template<class ForwardIterator>
inline size_t argmax(ForwardIterator first, ForwardIterator last)
{
    return std::distance(first, std::max_element(first, last));
}

GraphDeduplicator::GraphDeduplicator(vector<NODE> input){
    // 
    this->input = input;
    for(auto& node : input){
        for(auto& e : node.left){
            LM = max(LM, e + 1);
        }
        for(auto& e : node.right){
            RM = max(RM, e + 1);
        }
    }
    for (int i = 0; i < input.size(); ++i) {
        neighbour.push_back({});
    }
}

GraphDeduplicator::GraphDeduplicator(string filename) {
    // read file
    ifstream ifs("./trigraph/" + filename, std::ifstream::in);
    if (ifs) {
        string line;
        while (getline(ifs, line)) {
            ++N;
            if (line[0] == '#' || line[0] == '%') {
                continue;
            }
            vector<string> content = split_line(line, "\r\t\n ");
            int mid = atoi(content[0].c_str());
            int k = 1;
            int left_len = atoi(content[k++].c_str());
            set<int> left_tmp;
            for (int i = 0; i < left_len; i++) {
                int temp = atoi(content[k++].c_str());
                left_tmp.insert(temp);
                LM = max(LM, temp+1);
            }
            set<int> right_tmp;
            int right_len = atoi(content[k++].c_str());
            for (int i = 0; i < right_len; i++) {
                int temp = atoi(content[k++].c_str());
                right_tmp.insert(temp);
                RM = max(RM, temp+1);
            }
            // raw_edge_tot += (left_len + right_len);
            NODE node(mid, left_tmp, right_tmp);
            input.push_back(node);
        }
    }
    ifs.close();
    for (int i = 0; i < N; ++i) {
        neighbour.push_back({});
    }
    
}

inline set<int> GraphDeduplicator::findCommonNeighbour(const set<int>& a, const set<int>& b) {
    set<int> v;
    set_intersection(a.begin(), a.end(), b.begin(), b.end(), 
                    insert_iterator<set<int>>(v, v.begin()));
    return v;
}


void GraphDeduplicator::insert_edge(int u, int v, bool left){
    int v_offset = input.size() + v;
    if(!left){
        v_offset = input.size() + LM + v;
    }
    // assert(u < input.size() + LM + RM);
    // assert(v_offset < input.size() + LM + RM);

    next.push_back(Vnode[u]);
    if(left){
        edges.push_back(OUT_NODE(-1, v));
    }
    else{
        edges.push_back(OUT_NODE(1, v));
    }
    Vnode[u] = edges.size() - 1;

    next.push_back(Vnode[v_offset]);
    edges.push_back(OUT_NODE(0, u));
    Vnode[v_offset] = edges.size() - 1;
}


typedef unsigned int ui;

void GraphDeduplicator::deduplicateBySetCover(){
    printf("****************Set Cover****************\n");
    ui v_num = input.size() + RM + LM;
    Vnode = new ui[v_num];
    memset(Vnode, -1, sizeof(ui) * v_num);
    // construct the threepart graph
    for(ui s = input.size(), i = 0;i < s;i++){
        NODE& vn = input[i];
        for(auto& ln : vn.left){
            insert_edge(i, ln, 1);
        }
        for(auto& rn : vn.right){
            insert_edge(i, rn, 0);
        }
    }
    
    set<pair<pair<int, int>, int>> redundant_edges;
    // record weight for each edge.
    ui* counter = new ui[edges.size()];
    memset(counter, 0, sizeof(ui) * edges.size());

    ui* left_tag = new ui[LM];
    ui* right_tag = new ui[RM];
    memset(left_tag, 0, sizeof(ui) * LM);
    memset(right_tag, 0, sizeof(ui) * RM);

    for(int i = 0, s = input.size();i < s;i++){
        NODE& node_i = input[i];
        for(int j = i + 1;j < s;j++){
            NODE& node_j = input[j];
            if (*(node_j.left.begin()) > *(--node_i.left.end())) {
                break;
            }
            
            // detection for conflict can be improved.
            set<int> node_i_j_left = findCommonNeighbour(node_i.left, node_j.left);
            if (node_i_j_left.size() == 0) continue;
            set<int> node_i_j_right = findCommonNeighbour(node_i.right, node_j.right);
            if (node_i_j_right.size() == 0) continue;

            // count
            for(auto& l : node_i_j_left){
                left_tag[l] = 1;
            }
            for(auto& r : node_i_j_right){
                right_tag[r] = 1;
            }

            // traverse the adjacent list of vertex i and j.
            for(int e = Vnode[i]; e != -1; e = next[e]){
                if((edges[e].isVirtual == -1 && left_tag[edges[e].node_id] == 1 )){
                    counter[e] += node_i_j_right.size() ;
                } 
                else if(edges[e].isVirtual == 1 &&right_tag[edges[e].node_id] == 1 ){
                    counter[e] += node_i_j_left.size() ;
                }
            }
            
            for(int e = Vnode[j]; e != -1; e = next[e]){
                if((edges[e].isVirtual == -1 && left_tag[edges[e].node_id] == 1 )){
                    counter[e] += node_i_j_right.size() ;
                } 
                else if(edges[e].isVirtual == 1 &&right_tag[edges[e].node_id] == 1 ){
                    counter[e] += node_i_j_left.size() ;
                }
            }

            // recover the state
            for(auto& l : node_i_j_left){
                left_tag[l] = 0;
            }
            for(auto& r : node_i_j_right){
                right_tag[r] = 0;
            }
        }
    }

    ui conflicts = 0;
    for(int i = 0;i <edges.size();i++){
        conflicts += counter[i];
    }
    conflicts /= 4;
    printf("there are %d conflicts.\n", conflicts);
    
    ui* src = new ui[edges.size()];
    for(int i = 0, s = input.size() + RM + LM;i < s;i++){
        for(int e = Vnode[i];e != -1;e = next[e]){
            src[e] = i;
        }
    }

    ui max_freq = 0;
    for(int i = 0;i < edges.size();i++){
        max_freq = max(max_freq, counter[i]);
    }
    // delete the edge with highest score until no conflicts exist.
    
    ui* Order = new ui[max_freq + 1];
    memset(Order, -1, sizeof(ui) * (max_freq+1));
    ui* Onext = new ui[edges.size()];
    bool* IsEdge = new bool[edges.size()];
    memset(IsEdge, -1, sizeof(bool) * edges.size());

    for(int i = 0;i < input.size();i++){
        for(int e = Vnode[i]; e != -1; e = next[e]){
            Onext[e] = Order[counter[e]];
            Order[counter[e]] = e;
        }
    }
    
    int delete_edges = 0;
    ui times = 0;
    // resolve the conflicts until empty.
    for(int i = max_freq;i > 0;){
        while(i > 0 && Order[i] == -1){ i--;}
        if(i <= 0) break;
        int j = -1;
        // printf("i:%d \n", i);
        for(j = Order[i]; j != -1;){
            int tmp = Onext[j];
            if(IsEdge[j]) {
                if(counter[j] < i) {
                    Onext[j] = Order[counter[j]];
                    Order[counter[j]] = j;
                }
                else{
                    IsEdge[j] = 0;
                    // if(times <= 0) exit(1);
                    // printf("there are %d conflicts, %d left after this round.\n",conflicts, conflicts - counter[j]);
                    conflicts -= counter[j];
                    times--;
                    counter[j] = 0;
                    int s = src[j],  d = edges[j].node_id, type = edges[j].isVirtual;
                    ui* tag = nullptr;
                    delete_edges++;
                    if(type == -1){    // d on the left side
                        // printf("delete (%d, %d)L\n", s, d);
                        d += input.size();
                        tag = right_tag;
                    }
                    else if(type == 1){     //d on the right side
                        // printf("delete (%d, %d)R\n", s, d);
                        d += input.size() + LM;
                        tag = left_tag;
                    }

                    for(int e = Vnode[d]; e != -1;e = next[e]){
                        if(edges[e].isVirtual == 0 && edges[e].node_id == s){
                            IsEdge[e] = 0;
                            break;
                        }
                    }

                    for(int e = Vnode[s]; e != -1 ;e = next[e]){
                        if(edges[e].isVirtual == -type && IsEdge[e]){
                            tag[edges[e].node_id] = 1;
                        }
                    }
                    for(int e = Vnode[d]; e != -1;e = next[e]){
                        // check if this edge is valid.
                        // if it is virtual node.
                        if(IsEdge[e] && edges[e].isVirtual == 0 && edges[e].node_id != s){
                            int c = 0, t = -1;
                            for(int ve = Vnode[edges[e].node_id]; ve != -1;ve = next[ve]){
                                if(edges[ve].isVirtual == -type && IsEdge[ve] && tag[edges[ve].node_id] >= 1){
                                    counter[ve]--;
                                    tag[edges[ve].node_id]++;
                                    c++;
                                }
                                else if (edges[ve].isVirtual == type){
                                    if(type == -1 && (d - input.size() == edges[ve].node_id)){
                                        t = ve;
                                    }
                                    else if(type == 1 && (d - input.size() - LM == edges[ve].node_id)){
                                        t = ve;
                                    }
                                }
                            }
                            counter[t] -= c;
                        }
                    }
                    for(int e = Vnode[s]; e != -1;e = next[e]){
                        if(edges[e].isVirtual == -type && IsEdge[e]){
                            // add the broken connections.
                            if(tag[edges[e].node_id] == 1){
                                redundant_edges.insert(make_pair(make_pair(d, edges[e].node_id), type));
                            }
                            counter[e] = counter[e] - tag[edges[e].node_id] + 1;
                            tag[edges[e].node_id] = 0;
                        }
                    }
                    Order[i] = tmp;
                    break;
                }
            }
            j = tmp;
            // delete this edge and update the priority for each influenced edge.
        }
        if(j == -1) Order[i]= -1;
    }

    printf("delete edges:%d\n", delete_edges);
    // postprocess the degree one cases.
    for(int i = 0;i < input.size();i++){
        int left_nb = 0, right_nb = 0, l = 0, r = 0;
        for(int e = Vnode[i]; e != -1;e = next[e]){
            if(IsEdge[e]){
                if(edges[e].isVirtual == -1){
                    left_nb++;
                    l = edges[e].node_id;
                }
                else if(edges[e].isVirtual == 1){
                    right_nb++;
                    r = edges[e].node_id;
                }
            }
        }
        if(left_nb <= 1 ||  right_nb <= 1){
            // collapse the virtual i and add the direct edges.
            if(left_nb == 1){
                for(int e = Vnode[i]; e != -1;e = next[e]){
                    if(IsEdge[e] && edges[e].isVirtual == 1){
                        redundant_edges.insert(make_pair(make_pair(l + input.size(), edges[e].node_id), -1));
                    }
                }   
            }
            else if(right_nb == 1){
                for(int e = Vnode[i]; e != -1;e = next[e]){
                    if(IsEdge[e] && edges[e].isVirtual == -1){
                        redundant_edges.insert(make_pair(make_pair(r + input.size() + LM, edges[e].node_id), 1));
                    }
                }
            }
            for(int e = Vnode[i]; e != -1;e = next[e]){
                IsEdge[e] = 0;
            } 
        }
    }

    int c = 0;
    ui relations = 0;
    for(int i = 0;i < input.size();i++){
        int left_nb = 0, right_nb = 0;
        for(int e = Vnode[i]; e != -1;e = next[e]){
            if(IsEdge[e]){
                c++;
                if(edges[e].isVirtual == -1){
                    // printf("%d %dL\n",input[i].input_node_id, edges[e].node_id);
                    left_nb++;
                }
                else if(edges[e].isVirtual == 1){
                    // printf("%d %dR\n",input[i].input_node_id, edges[e].node_id);
                    right_nb++;
                }
            }
        }
        relations += left_nb * right_nb;
    }
    c += redundant_edges.size();
    relations += redundant_edges.size();
    
    printf("Number of edges after deduplication: %d\n", c);
    printf("Set cover Compression ratio: %f%\n", c * 100.0f / this->exp_size);
    printf("relations: %d\n", relations);

    delete[] IsEdge;
    delete[] src;
    delete[] Order;
    delete[] Onext;
    delete[] left_tag;
    delete[] right_tag;
    delete[] Vnode;
    delete[] counter;
}

void GraphDeduplicator::construct_graphTopology() {
    for (int i = 0; i < N; ++i) {
        NODE& node_i = input[i];
        for (int j = i + 1; j < N; ++j) {
            NODE& node_j = input[j];
            if (*(node_j.left.begin()) > *(--node_i.left.end())) {
                break;
            }
            // detection for conflict can be improved.
            set<int> node_i_j_left = findCommonNeighbour(node_i.left, node_j.left);
            if (node_i_j_left.size() == 0) continue;
            set<int> node_i_j_right = findCommonNeighbour(node_i.right, node_j.right);
            if (node_i_j_right.size() == 0) continue;
            neighbour[i].push_back(j);
            neighbour[j].push_back(i);
        }
    }
}

void GraphDeduplicator::assign_vertex_weight(int k) {

}


void GraphDeduplicator::build_conflict_graph() {
    // sort(input.begin(), input.end(), compare); 
    construct_graphTopology();
    weight.resize(N, 1);
    // assign_vertex_weight(1);
}

void GraphDeduplicator::dedup1(){
    build_revertedIDX_for_realNodes();
    
    set<pair<int,int>> redundant_edges;
    // set<pair<int,int>> redundant_edges;
    VVON tmp(out_left);
    // shuffle the vertices.
    srand(unsigned(time(0)));
    random_shuffle(tmp.begin(), tmp.end(), Rand);

    // random_shuffle();
    int s = 1;
    for(auto& real_node : tmp){
        // if(s == 1){
        //     // cout<< real_node << endl;
        //     // printf("%d \n", real_node[0].node_id);
        //     s++;
        // }
        set<int> processed;
        for(auto& vn : real_node){
            // some right neighbors of this virtual node are deleted in previous iterations.
            set<int> right_res = set_minus(input[vn.node_id].right, input[vn.node_id].right_del);
            set<int> intersec = set_intersec(right_res, processed);
            set<int> unique = set_minus(right_res, processed);
            // delete the common neighbors to avoid duplication.
            input[vn.node_id].right_del.insert(intersec.begin(), intersec.end());
            // add the missing direct edges between intersec and the left real nodes of this virtual node.
            for(auto& Ln : input[vn.node_id].left){
                for(auto& Rt : intersec){
                    redundant_edges.insert(make_pair(Ln, Rt));
                }
            }
            for(auto& Rn : unique){
                processed.insert(Rn);
            }
        }
        
    }

    out_left.clear();
    out_right.clear();
    for(auto& vn : input){
        int a = vn.left.size() - vn.left_del.size(), b = vn.right.size() - vn.right_del.size();
        if(a * b < a + b){
            set<int> left_res = set_minus(vn.left, vn.left_del), 
                    right_res = set_minus(vn.right, vn.right_del);
            vn.left_del.insert(vn.left.begin(), vn.left.end());
            vn.right_del.insert(vn.right.begin(), vn.right.end());
            for(auto& l : left_res){
                for(auto& r : right_res){
                    redundant_edges.insert(make_pair(l, r));
                }
            }
        }
    }

    build_revertedIDX_for_realNodes();

    // put the redundant edges back.
    for(auto& p : redundant_edges){
        int p1 = p.first, p2 = p.second;
        bool ok = 1;
        for(auto& n : out_left[p1]) if(n.isVirtual){
            set<int> right_res = set_minus(input[n.node_id].right, input[n.node_id].right_del);
            // form a triangle
            if(right_res.count(p2) > 0 ){
                ok = 0;
            }
        }
        if(ok){
            out_left[p1].push_back(OUT_NODE(0, p2));
            out_right[p2].push_back(OUT_NODE(0, p1));
        }
    }

    printf("****************DeDup1****************\n");
    report_result();

    for(auto& n : input){
        n.left_del.clear();
        n.right_del.clear();
    }
    out_left.clear();
    out_right.clear();
}

// greedy real node first deduplication
void GraphDeduplicator::greedyDedup(){
    printf("****************Greedy DeDup****************\n");
    build_revertedIDX_for_realNodes();
    
    // count the degree
    // clock_t start = clock();
    // int iterations = 1;
    // int count = 0;
    // int invalid = 0;
    // int total = 0;
    // while(iterations--){
    //     vector<int> flag(out_left.size(), 0);
    //     for(auto& l : out_left){
    //         vector<int> t;
    //         // set<int> tmp;
    //         for(auto& node : l){
    //             if(node.isVirtual){
    //                 for(auto& r : input[node.node_id].right){
    //                     // t.insert(r);
    //                     if(!flag[r]){
    //                         t.push_back(r);
    //                         flag[r] = 1;
    //                     }
    //                     // else{
    //                     //     invalid++;
    //                     // }
    //                     invalid++;
    //                     total++;
    //                 }
    //             }
    //             else{
    //                 if(!flag[node.node_id]){
    //                     t.push_back(node.node_id);
    //                     flag[node.node_id] = 1;
    //                 }
    //                 // else{
    //                 //     invalid++;
    //                 // }
    //                 // total++;
    //                 // t.push_back(node.node_id);
    //                 // count++;
    //             }
    //         }
    //         // 
    //         for(auto& x : t){
    //             flag[x] = 0;
    //         }
    //         total += t.size();
    //         // sort(t.begin(), t.end());
    //         // t.erase(unique(t.begin(), t.end()), t.end());
    //     }
    //     printf("invalid:%d total:%d ratio:%f\n",invalid, total, invalid * 1.0 / total);
    // }
    // clock_t end = clock();
    // printf("condensed graph count degree time cost:%fs count : %d \n", ((double)(end - start) / (CLOCKS_PER_SEC )), count);

    // k-core for duplicate graphs. 
    // start = clock();
    // int k = 4;
    // iterations =10;
    // while(iterations--){
    //     vector<int> k_list;
    //     vector<int> visited(out_left.size(), 0);
    //     vector<int> deg(out_left.size(), 0);
    //     vector<int> flag(out_left.size(), 0);   
    //     for(int i = 0, s = out_left.size();i < s;i++){
    //         auto& l = out_left[i]; 
            
    //         // count the degree.
    //         vector<int> t;
    //         // set<int> tmp;
    //         for(auto& node : l){
    //             if(node.isVirtual){
    //                 for(auto& r : input[node.node_id].right){
    //                     // t.insert(r);
    //                     // t.push_back(r);
    //                     // count++;
    //                     if(!flag[r]){
    //                         t.push_back(r);
    //                         flag[r] = 1;
    //                     }
    //                 }
    //             }
    //             else{
    //                 if(!flag[node.node_id]){
    //                     t.push_back(node.node_id);
    //                     flag[node.node_id] = 1;
    //                 }
                    
    //                 // t.push_back(node.node_id);
    //                 // count++;
    //             }
    //         }
    //         // sort(t.begin(), t.end());
    //         // t.erase(unique(t.begin(), t.end()), t.end());
    //         for(auto& x : t){
    //             flag[x] = 0;
    //         }
    //         deg[i] = t.size();
    //         if(count < k){
    //             k_list.push_back(i);
    //         }
    //     }   
        
    //     while(!k_list.empty()){
    //         int head = k_list.back();
    //         k_list.pop_back();
    //         // degree of head < k
    //         // remove this vertex and update its adjacent neighbors.
    //         visited[head] = 1;
    //         auto& l = out_left[head];    
    //         // count the degree.
            
    //         vector<int> t;
    //         // set<int> tmp;
    //         for(auto& node : l){
    //             if(node.isVirtual){
    //                 for(auto& r : input[node.node_id].right){
    //                     // t.insert(r);
    //                     // t.push_back(r);
    //                     // count++;
    //                     if(!flag[r]){
    //                         t.push_back(r);
    //                         flag[r] = 1;
    //                     }
    //                 }
    //             }
    //             else{
    //                 if(!flag[node.node_id]){
    //                     t.push_back(node.node_id);
    //                     flag[node.node_id] = 1;
    //                 }
    //                 // t.push_back(node.node_id);
    //                 // count++;
    //             }
    //         }
    //         // sort(t.begin(), t.end());
    //         // t.erase(unique(t.begin(), t.end()), t.end());

    //         for(auto& r : t){
    //             deg[r]--;
    //             if(!visited[r] && deg[r] < k){
    //                 k_list.push_back(r);
    //             }
    //         }
    //     }
    // }
    // end = clock();
    // printf("condensed k-core time cost:%fs \n", ((double)(end - start) / (CLOCKS_PER_SEC )));
    
    // bfs 
    // start = clock();
    // iterations =10;
    // while(iterations--){
    //     vector<int> k_list;
    //     vector<int> visited(out_left.size(), 0);
    //     // vector<int> deg(out_left.size(), 0);
    //     for(int i = 0, s = out_left.size();i < s;i++){
    //         if(!visited[i]){
    //             // 
    //             queue<int> q;
    //             q.push(i);
    //             visited[i] = 1;
    //             while(!q.empty()){
    //                 int v = q.front();
    //                 q.pop();
    //                 auto& l = out_left[v];    
    //                 // set<int> tmp;
    //                 for(auto& node : l){
    //                     if(node.isVirtual){
    //                         for(auto& r : input[node.node_id].right){
    //                             if(!visited[r]){
    //                                 visited[r] = 1;
    //                                 q.push(r);
    //                             }
    //                         }
    //                     }
    //                     else{
    //                         visited[node.node_id] = 1;
    //                         q.push(node.node_id);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    // end = clock();
    // printf("condensed bfs time cost:%fs \n", ((double)(end - start) / (CLOCKS_PER_SEC )));

    // pagerank
    // start = clock();
    // count = 0;
    // // printf("%d %d \n", LM, out_left.size());
    // vector<int> deg(LM, 0);
    // vector<int> flag(out_left.size(), 0);
    // for(auto& l : out_left){
    //     vector<int> t;
    //     // set<int> tmp;
    //     for(auto& node : l){
    //         if(node.isVirtual){
    //             for(auto& r : input[node.node_id].right){
    //                 // t.insert(r);
    //                 if(!flag[r]){
    //                     t.push_back(r);
    //                     flag[r] = 1;
    //                 }
                    
    //                 // count++;
    //             }
    //         }
    //         else{
    //             if(!flag[node.node_id]){
    //                 t.push_back(node.node_id);
    //                 flag[node.node_id] = 1;
    //             }
    //             // t.push_back(node.node_id);
    //             // count++;
    //         }
    //     }
    //     // 
    //     for(auto& x : t){
    //         flag[x] = 0;
    //     }
    //     deg[count++] = t.size();
    // }

    // // printf("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");

    // float* prev = new float[LM];
    // float* curr = new float[LM];
    // float damp = 0.85f;
    // iterations = 1;
    // assert(LM == out_left.size());
    // while(iterations--){
    //     for(int i = 0, s = out_left.size();i < s;i++){
    //         curr[i] = (1.0 - damp) / LM;
    //         float t = 0.0f;
    //         vector<int> visited(out_left.size(), 0);
    //         // printf("1\n");
    //         auto& l = out_left[i];    
    //         for(auto& node : l){
    //             if(node.isVirtual){
    //                 assert(node.node_id < input.size());
    //                 for(auto& r : input[node.node_id].right){
    //                     // printf("2\n");
    //                     assert(r < out_left.size());
    //                     assert(r >= 0);
    //                     if(!visited[r]){
    //                         assert(r < LM);
    //                         visited[r] = 1;
    //                         t +=  prev[r] / deg[r];
    //                     }
    //                 }
    //             }
    //             else{
    //                 assert(node.node_id < out_left.size());
    //                 assert(node.node_id >= 0);
    //                 if(!visited[node.node_id]){
    //                     assert(node.node_id < LM);
    //                     visited[node.node_id] = 1;
    //                     // q.push(node.node_id);
    //                     t += prev[node.node_id] / deg[node.node_id];
    //                 }
    //             }
    //             // printf("3\n");
    //         }
    //         curr[i] += damp * t;
    //     }
    //     swap(prev, curr);
    // }
    // end = clock();
    // printf("condensed pagerank time cost:%fs \n", ((double)(end - start) / (CLOCKS_PER_SEC)));

    set<pair<int,int>> redundant_edges;
    VVON tmp(out_left);
    // shuffle the vertices.
    srand(unsigned(time(0)));
    random_shuffle(tmp.begin(), tmp.end(), Rand);

    for(auto& real_node : tmp){
        set<int> processed;
        set<int> visited;
        int times = real_node.size();
        while(times--){
            int target = -1, max = -1;
            for(auto& vn : real_node){
                if(visited.count(vn.node_id) == 0){
                    set<int> right_res = set_minus(input[vn.node_id].right, input[vn.node_id].right_del);
                    // set<int> intersec = set_intersec(right_res, processed);
                    set<int> unique = set_minus(right_res, processed);
                    if((int)unique.size() > max){
                        max = unique.size();
                        target = vn.node_id;
                    }
                }
            }   
            visited.insert(target);
            set<int> right_res = set_minus(input[target].right, input[target].right_del);
            set<int> intersec = set_intersec(right_res, processed);
            set<int> unique = set_minus(right_res, processed);
            // delete the common neighbors to avoid duplication.
            input[target].right_del.insert(intersec.begin(), intersec.end());
            // add the missing direct edges between intersec and the left real nodes of this virtual node.
            for(auto& Ln : input[target].left){
                for(auto& Rt : intersec){
                    redundant_edges.insert(make_pair(Ln, Rt));
                }
            }
            for(auto& Rn : unique){
                processed.insert(Rn);
            }
        }
    }

    out_left.clear();
    out_right.clear();
    
    for(auto& vn : input){
        int a = vn.left.size() - vn.left_del.size(), b = vn.right.size() - vn.right_del.size();
        if(a * b < a + b){
            set<int> left_res = set_minus(vn.left, vn.left_del), 
                    right_res = set_minus(vn.right, vn.right_del);
            vn.left_del.insert(vn.left.begin(), vn.left.end());
            vn.right_del.insert(vn.right.begin(), vn.right.end());
            for(auto& l : left_res){
                for(auto& r : right_res){
                    redundant_edges.insert(make_pair(l, r));
                }
            }
        }
    }

    
    build_revertedIDX_for_realNodes();

    // put the redundant edges back.
    for(auto& p : redundant_edges){
        int p1 = p.first, p2 = p.second;
        bool ok = 1;
        for(auto& n : out_left[p1]) if(n.isVirtual){
            set<int> right_res = set_minus(input[n.node_id].right, input[n.node_id].right_del);
            // form a triangle
            if(right_res.count(p2) > 0 ){
                ok = 0;
            }
        }
        if(ok){
            out_left[p1].push_back(OUT_NODE(0, p2));
            out_right[p2].push_back(OUT_NODE(0, p1));
        }
    }

    report_result();

    // 因为在进行countdegree的时候不想每次重新求一次差集，所以直接求一次。但之后要记得复原。
    vector<NODE> input_copy = input;
    for(auto& vn : input){
        set<int> left_res = set_minus(vn.left, vn.left_del), 
                right_res = set_minus(vn.right, vn.right_del);
        vn.left = left_res;
        vn.right = right_res;
    }   
    
    // count the degree.
    // start = clock();
    // iterations =10;
    // count = 0;
    // while(iterations--){
    //     for(auto& l : out_left){
    //         for(auto& node : l){
    //             if(node.isVirtual){
    //                 for(auto& r : input[node.node_id].right){
    //                     count++;
    //                 }
    //             }
    //             else{
    //                 count++;
    //             }
    //         }
    //     }
    // }
    // end = clock();
    // printf("greedy Deduplication s degree time cost:%fs  count : %d \n", ((double)(end - start) / (CLOCKS_PER_SEC )), count);

    // // k-core for non-duplicate graphs.
    // start = clock();
    // k = 4;
    // iterations =10;
    // while(iterations--){
    //     vector<int> k_list;
    //     vector<int> visited(out_left.size(), 0);
    //     vector<int> deg(out_left.size(), 0);
    //     for(int i = 0, s = out_left.size();i < s;i++){
    //         auto& l = out_left[i];    
    //         // count the degree.
    //         for(auto& node : l){
    //             if(node.isVirtual){
    //                 for(auto& r : input[node.node_id].right){
    //                     count++;
    //                 }
    //             }
    //             else{
    //                 count++;
    //             }
    //         }
    //         deg[i] = count;
    //         if(count < k){
    //             k_list.push_back(i);
    //         }
    //     }   
        
    //     while(!k_list.empty()){
    //         int head = k_list.back();
    //         k_list.pop_back();
    //         // degree of head < k
    //         // remove this vertex and update its adjacent neighbors.
    //         visited[head] = 1;
    //         auto& l = out_left[head];    
    //         // count the degree.
            
    //         for(auto& node : l){
    //             if(node.isVirtual){
    //                 for(auto& r : input[node.node_id].right){
    //                     // count++;
    //                     deg[r]--;
    //                     if(!visited[r] && deg[r] < k){
    //                         k_list.push_back(r);
    //                     }
    //                 }
    //             }
    //             else{
    //                 int r = node.node_id;
    //                 deg[r]--;
    //                 if(!visited[r] && deg[r] < k){
    //                     k_list.push_back(r);
    //                 }
    //             }
    //         }
    //     }
    // }
    // end = clock();
    // printf("greedy Deduplications k-core time cost:%fs \n", ((double)(end - start) / (CLOCKS_PER_SEC )));

    // // bfs
    // start = clock();
    // iterations = 1;
    // while(iterations--){
    //     vector<int> visited(out_left.size(), 0);
    //     // vector<int> deg(out_left.size(), 0);
    //     for(int i = 0, s = out_left.size();i < s;i++){
    //         if(!visited[i]){
    //             // 
    //             queue<int> q;
    //             q.push(i);
    //             visited[i] = 1;
    //             while(!q.empty()){
    //                 int v = q.front();
    //                 q.pop();
    //                 // printf("a\n");
    //                 auto& l = out_left[v];    
    //                 // set<int> tmp;
    //                 for(auto& node : l){
    //                     if(node.isVirtual){
    //                         for(auto& r : input[node.node_id].right){
    //                             if(!visited[r]){
    //                                 visited[r] = 1;
    //                                 q.push(r);
    //                             }
    //                         }
    //                     }
    //                     else{
    //                         if(!visited[node.node_id]){
    //                             visited[node.node_id] = 1;
    //                             q.push(node.node_id);
    //                         }
                            
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    // end = clock();
    // printf("greedy deduplication bfs time cost:%fs \n", ((double)(end - start) / (CLOCKS_PER_SEC )));

    // // pagerank
    // start = clock();
    // count = 0;
    // // vector<int> deg(LM, 0);
    
    // for(auto& l : out_left){
    //     vector<int> t;
    //     // set<int> tmp;
    //     for(auto& node : l){
    //         if(node.isVirtual){
    //             for(auto& r : input[node.node_id].right){
    //                 // t.insert(r);
    //                 t.push_back(r);
    //                 flag[r] = 1;
    //             }
    //         }
    //         else{
    //             t.push_back(node.node_id);
    //             // t.push_back(node.node_id);
    //             // count++;
    //         }
    //     }
    //     deg[count++] = t.size();
    // }

    // // float* prev = new float[LM];
    // // float* curr = new float[LM];
    // damp = 0.85f;
    // iterations = 1;
    // while(iterations--){
    //     for(int i = 0, s = out_left.size();i < s;i++){
    //         curr[i] = (1 - damp) / LM;
    //         float t = 0.0f;
    //         auto& l = out_left[i];    
    //         for(auto& node : l){
    //             if(node.isVirtual){
    //                 for(auto& r : input[node.node_id].right){
    //                     t +=  prev[r] / deg[r];
    //                 }
    //             }
    //             else{
    //                 t += prev[node.node_id] / deg[node.node_id];
    //             }
    //         }
    //     }
    //     swap(prev, curr);
    // }
    // end = clock();
    // printf("condensed pagerank time cost:%fs \n", ((double)(end - start) / (CLOCKS_PER_SEC )));

    // delete[] prev;
    // delete[] curr;
    input = input_copy;
    // clean up 
    for(auto& n : input){
        n.left_del.clear();
        n.right_del.clear();
    }
    out_left.clear();
    out_right.clear();
}


void GraphDeduplicator::bfsTest(int iter){
    
}


void GraphDeduplicator::count_expand_edges(){
    vector<vector<int>> tmp_left;
    tmp_left.resize(LM + 1, {});
    for(int s = input.size(), i = 0;i < s;i++){
        NODE& node = input[i];
        for(auto& Lrn : node.left){
            assert(Lrn <= LM);
            tmp_left[Lrn].push_back(i);
        }
    }

    unsigned int counter = 0;
    for(auto& x : tmp_left){
        set<int> visited;
        for(auto& mid : x){
            for(auto& r : input[mid].right){
                visited.insert(r);
            }
        }
        counter += visited.size();
    }

    this->exp_size = counter;
    printf("size of expanded graph: %d \n", counter);
}



void GraphDeduplicator::remove_edge_and_restore_connections(int target, bool left, set<int>& removed_edges, set<int>& common, set<pair<int,int>>& redundant_edges){
    // left side
    if(left){
        input[target].left_del.insert(removed_edges.begin(), removed_edges.end());
        // add the connections between left common and right unique.
        set<int> unique = set_minus(input[target].right, common);
        for(auto& i : removed_edges) {
            for(auto& j : unique) {
                redundant_edges.insert(make_pair(i, j));
            }    
        }
    }
    // right side
    else{
        input[target].right_del.insert(removed_edges.begin(), removed_edges.end());
        // add the connections between right common and left unique.
        set<int> unique = set_minus(input[target].left, common);
        for(auto& i : removed_edges) {
            for(auto& j : unique) {
                redundant_edges.insert(make_pair(j, i));
            }    
        }
    }
}

void GraphDeduplicator::build_revertedIDX_for_realNodes(){
    for (int i = 0; i < LM; ++i) {
        out_left.push_back({});
    }
    for (int i = 0; i < RM; ++i) {
        out_right.push_back({});
    }
    for(int s = input.size(), i = 0;i < s;i++){
        NODE& node = input[i];
        set<int> left_residual = set_minus(node.left, node.left_del), right_residual = set_minus(node.right, node.right_del);
        if(left_residual.size() > 0 && right_residual.size() > 0){
            for(auto& Lrn : left_residual){
                out_left[Lrn].push_back(OUT_NODE(1, i));
            }
            for(auto & Rrn : right_residual){
                out_right[Rrn].push_back(OUT_NODE(1, i));
            }
        }
    }
}


void GraphDeduplicator::deduplicateBySearch(const vector<int>& mvc, const vector<int>& mis) {
    printf("****************Search****************\n");
    vector<int> unconflict_set;
    unconflict_set.resize(input.size(), 0);
    set<pair<int,int>> redundant_edges;

    // label each vertex which is mutually unconflict.
    for(auto n : mis){
        unconflict_set[n] = 1;
    }

    int valid_check = 0;
    long long int cost = 0;
    // deduplicate one by one
    for(auto id : mvc){
        for(auto& ngb : neighbour[id]){
            // conflict exist, resolve it
            if(unconflict_set[ngb] == 1){
                valid_check++;
                // find left and right common neighbors
                set<int> left_residual = set_minus(input[ngb].left,  input[ngb].left_del), 
                         right_residual = set_minus(input[ngb].right, input[ngb].right_del);
                set<int> left_common = findCommonNeighbour(input[id].left, left_residual);
                set<int> right_common = findCommonNeighbour(input[id].right, right_residual);
                
                // choose which side to handle 
                int id_unique_left = input[id].left.size() - left_common.size();
                int id_unique_right = input[id].right.size() - right_common.size();
                int ngb_unique_left = input[ngb].left.size() - left_common.size();
                int ngb_unique_right = input[ngb].right.size() - right_common.size();
                
                vector<int> choice = { (int)left_common.size() * (1 - id_unique_right), 
                                    (int)left_common.size() * (1 - ngb_unique_right),
                                    (int)right_common.size() * (1 - id_unique_left),
                                    (int)right_common.size() * (1 - ngb_unique_left)};
                int idx = argmax(choice.begin(), choice.end());

                cost += choice[idx];
                // vector<int> choice = { (int)left_common.size() * id_unique_right, 
                //                     (int)left_common.size() * ngb_unique_right,
                //                     (int)right_common.size() * id_unique_left,
                //                     (int)right_common.size() * ngb_unique_left};
                // int idx = argmin(choice.begin(), choice.end());
                int target_node = id;
                set<int> removed_edge, oppo_common;
                bool left = 1;
                // arguments initialization
                switch(idx){
                    case 0: 
                        target_node = id, removed_edge = left_common, oppo_common = right_common; 
                        break;
                    case 1: 
                        target_node = ngb, removed_edge = left_common, oppo_common = right_common; 
                        break;
                    case 2: 
                        target_node = id, removed_edge = right_common, oppo_common = left_common, left = 0; 
                        break;
                    case 3: 
                        target_node = ngb, removed_edge = right_common, oppo_common = left_common, left = 0; 
                        break;
                }
                // remove the edge to deduplicate
                // restore the destroyed connections caused by deduplications.
                remove_edge_and_restore_connections(target_node, left, removed_edge, oppo_common, redundant_edges);
            }
        }
        unconflict_set[id] = 1;
    }

    printf("valid check: %d, predictive cost:%lld\n", valid_check, cost);
    // postprocess the virtual vertices.
    // for virtual vertices where ab < a+b, we remove it.
    for(auto& vn : input){
        int a = vn.left.size() - vn.left_del.size(), b = vn.right.size() - vn.right_del.size();
        if(a * b < a + b){
            set<int> left_res = set_minus(vn.left, vn.left_del), 
                    right_res = set_minus(vn.right, vn.right_del);
            vn.left_del.insert(vn.left.begin(), vn.left.end());
            vn.right_del.insert(vn.right.begin(), vn.right.end());
            for(auto& l : left_res){
                for(auto& r : right_res){
                    redundant_edges.insert(make_pair(l, r));
                }
            }
        }
    }

    build_revertedIDX_for_realNodes();

    // NOTE: need to construct the reverted index for real nodes.
    // to check if there is triangles. if no, the edge should be preserved, if any, it shoulde be removed.

    char* tag = new char[input.size()];
    memset(tag, 0, sizeof(char) * input.size());
    for(auto& p : redundant_edges){
        int p1 = p.first, p2 = p.second;
        for(auto n : out_left[p1]) if(n.isVirtual){
            tag[n.node_id] = 1;
        }
        bool ok = 1;
        for(auto n :  out_right[p2]) if(n.isVirtual){
            if(tag[n.node_id]){
                ok = 0;
            }
        }
        if(ok){
            out_left[p1].push_back(OUT_NODE(0, p2));
            out_right[p2].push_back(OUT_NODE(0, p1));
        }
        for(auto n : out_left[p1]) if(n.isVirtual){
            tag[n.node_id] = 0;
        } 
    }

    report_result();

    // 因为在进行countdegree的时候不想每次重新求一次差集，所以直接求一次。但之后要记得复原。
    vector<NODE> input_copy = input;
    for(auto& vn : input){
        set<int> left_res = set_minus(vn.left, vn.left_del), 
                right_res = set_minus(vn.right, vn.right_del);
        vn.left = left_res;
        vn.right = right_res;
    }   
    vector<int> deg(LM, 0);
    vector<int> flag(out_left.size(), 0);
    
    // count the degree.
    // clock_t start = clock();
    // int iterations =10;
    // int count = 0;
    // while(iterations--){
    //     for(auto& l : out_left){
    //         for(auto& node : l){
    //             if(node.isVirtual){
    //                 for(auto& r : input[node.node_id].right){
    //                     count++;
    //                 }
    //             }
    //             else{
    //                 count++;
    //             }
    //         }
    //     }
    // }
    // clock_t end = clock();
    // printf("greedy Deduplication s degree time cost:%fs  count : %d \n", ((double)(end - start) / (CLOCKS_PER_SEC )), count);

    // // k-core for non-duplicate graphs.
    // start = clock();
    // int k = 4;
    // iterations =10;
    // while(iterations--){
    //     vector<int> k_list;
    //     vector<int> visited(out_left.size(), 0);
    //     vector<int> deg(out_left.size(), 0);
    //     for(int i = 0, s = out_left.size();i < s;i++){
    //         auto& l = out_left[i];    
    //         // count the degree.
    //         for(auto& node : l){
    //             if(node.isVirtual){
    //                 for(auto& r : input[node.node_id].right){
    //                     count++;
    //                 }
    //             }
    //             else{
    //                 count++;
    //             }
    //         }
    //         deg[i] = count;
    //         if(count < k){
    //             k_list.push_back(i);
    //         }
    //     }   
        
    //     while(!k_list.empty()){
    //         int head = k_list.back();
    //         k_list.pop_back();
    //         // degree of head < k
    //         // remove this vertex and update its adjacent neighbors.
    //         visited[head] = 1;
    //         auto& l = out_left[head];    
    //         // count the degree.
            
    //         for(auto& node : l){
    //             if(node.isVirtual){
    //                 for(auto& r : input[node.node_id].right){
    //                     // count++;
    //                     deg[r]--;
    //                     if(!visited[r] && deg[r] < k){
    //                         k_list.push_back(r);
    //                     }
    //                 }
    //             }
    //             else{
    //                 int r = node.node_id;
    //                 deg[r]--;
    //                 if(!visited[r] && deg[r] < k){
    //                     k_list.push_back(r);
    //                 }
    //             }
    //         }
    //     }
    // }
    // end = clock();
    // printf("greedy Deduplications k-core time cost:%fs \n", ((double)(end - start) / (CLOCKS_PER_SEC )));

    // // bfs
    // start = clock();
    // iterations =10;
    // while(iterations--){
    //     vector<int> visited(out_left.size(), 0);
    //     // vector<int> deg(out_left.size(), 0);
    //     for(int i = 0, s = out_left.size();i < s;i++){
    //         if(!visited[i]){
    //             // 
    //             queue<int> q;
    //             q.push(i);
    //             visited[i] = 1;
    //             while(!q.empty()){
    //                 int v = q.front();
    //                 q.pop();
    //                 // printf("a\n");
    //                 auto& l = out_left[v];    
    //                 // set<int> tmp;
    //                 for(auto& node : l){
    //                     if(node.isVirtual){
    //                         for(auto& r : input[node.node_id].right){
    //                             if(!visited[r]){
    //                                 visited[r] = 1;
    //                                 q.push(r);
    //                             }
    //                         }
    //                     }
    //                     else{
    //                         if(!visited[node.node_id]){
    //                             visited[node.node_id] = 1;
    //                             q.push(node.node_id);
    //                         }
                            
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    // end = clock();
    // printf("greedy deduplication bfs time cost:%fs \n", ((double)(end - start) / (CLOCKS_PER_SEC )));

    // // pagerank
    // start = clock();
    // count = 0;
    // // vector<int> deg(LM, 0);
    
    // for(auto& l : out_left){
    //     vector<int> t;
    //     // set<int> tmp;
    //     for(auto& node : l){
    //         if(node.isVirtual){
    //             for(auto& r : input[node.node_id].right){
    //                 // t.insert(r);
    //                 t.push_back(r);
    //                 flag[r] = 1;
    //             }
    //         }
    //         else{
    //             t.push_back(node.node_id);
    //             // t.push_back(node.node_id);
    //             // count++;
    //         }
    //     }
    //     deg[count++] = t.size();
    // }

    // float* prev = new float[LM];
    // float* curr = new float[LM];
    // float damp = 0.85f;
    // iterations = 1;
    // while(iterations--){
    //     for(int i = 0, s = out_left.size();i < s;i++){
    //         curr[i] = (1 - damp) / LM;
    //         float t = 0.0f;
    //         auto& l = out_left[i];    
    //         for(auto& node : l){
    //             if(node.isVirtual){
    //                 for(auto& r : input[node.node_id].right){
    //                     t +=  prev[r] / deg[r];
    //                 }
    //             }
    //             else{
    //                 t += prev[node.node_id] / deg[node.node_id];
    //             }
    //         }
    //     }
    //     swap(prev, curr);
    // }
    // end = clock();
    // printf("condensed pagerank time cost:%fs \n", ((double)(end - start) / (CLOCKS_PER_SEC )));
    // delete[] tag;
    // delete[] prev;
    // delete[] curr;

    input = input_copy;

    // print_graph();
}

// advanced search, combine set cover and search method.
void GraphDeduplicator::deduplicateByWeightedSetCover(){
    printf("****************Greeedy Weightd Set Search****************\n");

    ui v_num = input.size() + RM + LM;
    Vnode = new ui[v_num];
    memset(Vnode, -1, sizeof(ui) * v_num);
    // construct the threepart graph
    for(ui s = input.size(), i = 0;i < s;i++){
        NODE& vn = input[i];
        for(auto& ln : vn.left){
            insert_edge(i, ln, 1);
        }
        for(auto& rn : vn.right){
            insert_edge(i, rn, 0);
        }
    }
    
    set<pair<pair<int, int>, int>> redundant_edges;
    // record weight for each edge.
    ui* counter = new ui[edges.size()];
    memset(counter, 0, sizeof(ui) * edges.size());

    ui* cost = new ui[edges.size()];
    memset(cost, 0, sizeof(ui) * edges.size());

    ui* left_tag = new ui[LM];
    ui* right_tag = new ui[RM];
    memset(left_tag, 0, sizeof(ui) * LM);
    memset(right_tag, 0, sizeof(ui) * RM);

    for(int i = 0, s = input.size();i < s;i++){
        NODE& node_i = input[i];
        for(int j = i + 1;j < s;j++){
            NODE& node_j = input[j];
            if (*(node_j.left.begin()) > *(--node_i.left.end())) {
                break;
            }
            
            // detection for conflict can be improved.
            set<int> node_i_j_left = findCommonNeighbour(node_i.left, node_j.left);
            if (node_i_j_left.size() == 0) continue;
            set<int> node_i_j_right = findCommonNeighbour(node_i.right, node_j.right);
            if (node_i_j_right.size() == 0) continue;

            // count
            for(auto& l : node_i_j_left){
                left_tag[l] = 1;
            }
            for(auto& r : node_i_j_right){
                right_tag[r] = 1;
            }

            // traverse the adjacent list of vertex i and j.
            for(int e = Vnode[i]; e != -1; e = next[e]){
                if((edges[e].isVirtual == -1 && left_tag[edges[e].node_id] == 1 )){
                    counter[e] += node_i_j_right.size() ;
                } 
                else if(edges[e].isVirtual == 1 && right_tag[edges[e].node_id] == 1 ){
                    counter[e] += node_i_j_left.size() ;
                }
            }
            
            for(int e = Vnode[j]; e != -1; e = next[e]){
                if((edges[e].isVirtual == -1 && left_tag[edges[e].node_id] == 1 )){
                    counter[e] += node_i_j_right.size() ;
                } 
                else if(edges[e].isVirtual == 1 && right_tag[edges[e].node_id] == 1 ){
                    counter[e] += node_i_j_left.size() ;
                }
            }

            // recover the state
            for(auto& l : node_i_j_left){
                left_tag[l] = 0;
            }
            for(auto& r : node_i_j_right){
                right_tag[r] = 0;
            }
        }
    }

    #ifdef DEBUG
    ui conflicts = 0;
    for(int i = 0;i < edges.size();i++){
        conflicts += counter[i];
    }
    assert(conflicts % 4 == 0);
    conflicts /= 4;
    printf("there are %d conflicts.\n", conflicts);
    #endif

    ui* src = new ui[edges.size()];
    for(int i = 0, s = input.size() + RM + LM;i < s;i++){
        for(int e = Vnode[i];e != -1;e = next[e]){
            src[e] = i;
        }
    }

    for(int i = input.size(), k = input.size() + LM ;i < k;i++){
        for(int e = Vnode[i];e != -1;e = next[e]){
            ui mid = edges[e].node_id;
            // tag the vertices visited on right side.
            for(int ve = Vnode[mid]; ve != -1;ve = next[ve]){
                if(edges[ve].isVirtual == 1){
                    right_tag[edges[ve].node_id] += 1;
                }   
            }
        }
        // tag the vertices visited on right side.
        
        for(int e = Vnode[i];e != -1;e = next[e]){
            ui mid = edges[e].node_id;
            ui c = 0;
            int t = -1;
            for(int ve = Vnode[mid]; ve != -1;ve = next[ve]){
                if(edges[ve].isVirtual == 1 && right_tag[edges[ve].node_id] == 1){
                    if(counter[ve] > 0){
                        cost[ve]++;
                    }
                    c++;
                } 
                else if(edges[ve].isVirtual == -1 && edges[ve].node_id == i - input.size()){
                    t = ve;
                }
            }   
            if(counter[t] > 0){
                cost[t] = c;
            }
        }
        for(int e = Vnode[i];e != -1;e = next[e]){
            ui mid = edges[e].node_id;
            // tag the vertices visited on right side.
            for(int ve = Vnode[mid]; ve != -1;ve = next[ve]){
                if(edges[ve].isVirtual == 1){
                    right_tag[edges[ve].node_id] = 0;
                }   
            }
        }
    }

    priority_queue<pair<float, int>> pq;
    ui* tag = nullptr;
    for(int i = 0;i < edges.size();i++){
        if(counter[i] > 0){
            // printf("%d %d %d, %d\n", src[i], edges[i].node_id, edges[i].isVirtual, cost[i]);
            pq.push(make_pair( cost[i] * -1.0 / counter[i] , i));
        }
    }

    ui max_freq = 0;
    for(int i = 0;i < edges.size();i++){
        max_freq = max(max_freq, counter[i]);
    }
    // delete the edge with highest score until no conflicts exist.


    bool* IsEdge = new bool[edges.size()];
    memset(IsEdge, -1, sizeof(bool) * edges.size());

    int delete_edges = 0;
    // resolve the conflicts until empty.
    while(!pq.empty()){
        if(conflicts == 0) {
            break;
        } 
        pair<float, int> p = pq.top();
        // printf("%f \n", p.first);
        pq.pop();
        int j = p.second;
        // pop the edge with the highest score and update others' score.
        if(IsEdge[j] != 0 && counter[j] > 0 && fabsf(p.first - cost[j] * -1.0f / counter[j]) < EPSINON) {
            // printf("%f %f\n", p.first,  (cost[j] * -1.0 / counter[j]));
            IsEdge[j] = 0;
            conflicts -= counter[j];
            counter[j] = 0;
            int s = src[j],  d = edges[j].node_id, type = edges[j].isVirtual;
            ui* tag = nullptr;
            delete_edges++;
            if(type == -1){    // d on the left side
                // printf("delete (%d, %d)L\n", s, d);
                d += input.size();
                tag = right_tag;
            }
            else if(type == 1){     //d on the right side
                // printf("delete (%d, %d)R\n", s, d);
                d += input.size() + LM;
                tag = left_tag;
            }

            for(int e = Vnode[d]; e != -1;e = next[e]){
                if(edges[e].isVirtual == 0 && edges[e].node_id == s){
                    IsEdge[e] = 0;
                    break;
                }
            }

            for(int e = Vnode[s]; e != -1 ;e = next[e]){
                if(edges[e].isVirtual == -type && IsEdge[e]){
                    tag[edges[e].node_id] = 1;
                }
            }
            for(int e = Vnode[d]; e != -1;e = next[e]){
                // check if this edge is valid.
                // if it is virtual node.
                if(IsEdge[e] && edges[e].isVirtual == 0 && edges[e].node_id != s){
                    int c = 0, t = -1;
                    for(int ve = Vnode[edges[e].node_id]; ve != -1;ve = next[ve]){
                        if(edges[ve].isVirtual == -type && IsEdge[ve] && tag[edges[ve].node_id] >= 1){
                            counter[ve]--;
                            tag[edges[ve].node_id]++;
                            c++;
                        }
                        else if (edges[ve].isVirtual == type){
                            if(type == -1 && (d - input.size() == edges[ve].node_id)){
                                t = ve;
                            }
                            else if(type == 1 && (d - input.size() - LM == edges[ve].node_id)){
                                t = ve;
                            }
                        }
                    }
                    counter[t] -= c;
                }
            }
            for(int e = Vnode[d]; e != -1;e = next[e]){
                // if it is virtual node.
                if(IsEdge[e] && edges[e].isVirtual == 0 && edges[e].node_id != s){
                    int c = 0, t = -1;
                    for(int ve = Vnode[edges[e].node_id]; ve != -1;ve = next[ve]){
                        if(edges[ve].isVirtual == -type && IsEdge[ve] && tag[edges[ve].node_id] >= 1){
                            if(tag[edges[ve].node_id] == 2) {
                                cost[ve]++;  
                                c++; 
                            }
                            if(counter[ve] > 0){
                                pq.push(make_pair(cost[ve] * -1.0f / counter[ve] , ve));
                            }
                        }
                        else if (edges[ve].isVirtual == type){
                            if(type == -1 && (d - input.size() == edges[ve].node_id)){
                                t = ve;
                            }
                            else if(type == 1 && (d - input.size() - LM == edges[ve].node_id)){
                                t = ve;
                            }
                        }
                    }
                    cost[t] += c;
                    if(counter[t] > 0){
                        pq.push(make_pair(cost[t] * -1.0f / counter[t], t));
                    }
                }
            }
            for(int e = Vnode[s]; e != -1;e = next[e]){
                if(edges[e].isVirtual == -type && IsEdge[e]){
                    // add the broken connections.
                    counter[e] = counter[e] - tag[edges[e].node_id] + 1;
                    if(tag[edges[e].node_id] == 1){
                        cost[e]--;
                        redundant_edges.insert(make_pair(make_pair(d, edges[e].node_id), type));
                    }
                    if(counter[e] > 0){
                        pq.push(make_pair(cost[e] * -1.0f / counter[e], e));
                    }
                    tag[edges[e].node_id] = 0;
                }
            }
        }
    }

    for(int i = 0;i < edges.size();i++){
        if(counter[i] > 0){
            printf("error\n");
        }
    }

    printf("delete edges:%d\n", delete_edges);

    // postprocess the degree one cases.
    for(int i = 0;i < input.size();i++){
        int left_nb = 0, right_nb = 0, l = 0, r = 0;
        for(int e = Vnode[i]; e != -1;e = next[e]){
            if(IsEdge[e]){
                if(edges[e].isVirtual == -1){
                    left_nb++;
                    l = edges[e].node_id;
                }
                else if(edges[e].isVirtual == 1){
                    right_nb++;
                    r = edges[e].node_id;
                }
            }
        }
        if(left_nb <= 1 ||  right_nb <= 1){
            // collapse the virtual i and add the direct edges.
            if(left_nb == 1){
                for(int e = Vnode[i]; e != -1;e = next[e]){
                    if(IsEdge[e] && edges[e].isVirtual == 1){
                        redundant_edges.insert(make_pair(make_pair(l + input.size(), edges[e].node_id), -1));
                    }
                }   
            }
            else if(right_nb == 1){
                for(int e = Vnode[i]; e != -1;e = next[e]){
                    if(IsEdge[e] && edges[e].isVirtual == -1){
                        redundant_edges.insert(make_pair(make_pair(r + input.size() + LM, edges[e].node_id), 1));
                    }
                }
            }
            for(int e = Vnode[i]; e != -1;e = next[e]){
                IsEdge[e] = 0;
            } 
        }
    }

    int c = 0;
    ui relations = 0;
    for(int i = 0;i < input.size();i++){
        int left_nb = 0, right_nb = 0;
        for(int e = Vnode[i]; e != -1;e = next[e]){
            if(IsEdge[e]){
                c++;
                if(edges[e].isVirtual == -1){
                    left_nb++;
                }
                else if(edges[e].isVirtual == 1){
                    right_nb++;
                }
            }
        }
        relations += left_nb * right_nb;
    }
    c += redundant_edges.size();
    relations += redundant_edges.size();

    printf("Number of edges after deduplication: %d\n", c);
    printf("Greedy Weighted Set cover Compression ratio: %f%\n", c * 100.0 / this->exp_size);
    printf("relations: %d\n", relations);

    // build list for direct edges.
    vector<int> left[LM], right[RM];
    for(auto& pp : redundant_edges){
        int type = pp.second;
        int u = pp.first.first, v = pp.first.second;
        if(type == -1){
            u -= input.size();
            assert(u < LM);
            left[u].push_back(v);
        }
        else if(type == 1){
            u -= (input.size() + LM);
            assert(u < RM);
            right[u].push_back(v);
        }
    }

    
    vector<int> left_neighbor_size(input.size() + 1, 0);
    vector<int> right_neighbor_size(input.size() + 1, 0);
    
    // count the neighbor size of each mid node.
    for(int i = 0;i < input.size();i++){
        int l = 0, r = 0;
        for(int e = Vnode[i]; e != -1;e = next[e]){
            if(IsEdge[e]){
                int type = edges[e].isVirtual;
                if(type == -1){
                    l++;
                }
                else{
                    r++;
                }
            }
        }
        left_neighbor_size[i] = l;
        right_neighbor_size[i] = r;
    }

    // // count degree
    // clock_t start = clock();
    // int iterations =10;
    // int count = 0;
    // while(iterations--){
    //     for(int i = input.size(), s = input.size() + LM;i < s;i++){
    //         for(int e = Vnode[i];e != -1; e = next[e]){
    //             if(IsEdge[e]){
    //                 for(int ee = Vnode[edges[e].node_id]; ee!= -1; ee = next[ee]){
    //                     if(IsEdge[ee] && edges[ee].isVirtual == 1){
    //                         count++;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     for(auto& p: redundant_edges){
    //        count++;
    //     }
    // }
    // clock_t end = clock();
    // printf("Weighted set cover count degree time:%f s count : %d \n", ((double)(end - start) / (CLOCKS_PER_SEC )), count);

    // // k-core
    // start = clock();
    // int k = 4;
    // iterations =10;
    // while(iterations--){
    //     vector<int> deg(LM, 0);
    //     vector<int> k_list;
    //     vector<int> visited(LM, 0);
    //     for(int i = input.size(), s = input.size() + LM;i < s;i++){
    //         int count = 0;
    //         for(int e = Vnode[i];e != -1; e = next[e]){
    //             if(IsEdge[e]){
    //                 for(int ee = Vnode[edges[e].node_id]; ee!= -1; ee = next[ee]){
    //                     if(IsEdge[ee] && edges[ee].isVirtual == 1){
    //                         count++;
    //                     }
    //                 }
    //             }
    //         }
    //         if(left[i-input.size()].size() + count < k){
    //             k_list.push_back(i-input.size());
    //             visited[i-input.size()] = 1;
    //         }
    //         deg[i-input.size()] = left[i-input.size()].size() + count;
    //     }

    //     while(!k_list.empty()){
    //         int head = k_list.back();
    //         k_list.pop_back();
    //         // visited[head] = 1;
    //         for(int e = Vnode[head];e != -1; e = next[e]){
    //             if(IsEdge[e]){
    //                 for(int ee = Vnode[edges[e].node_id]; ee!= -1; ee = next[ee]){
    //                     if(IsEdge[ee] && edges[ee].isVirtual == 1 && !visited[edges[ee].node_id]){
    //                         deg[edges[ee].node_id]--;
    //                         if(deg[edges[ee].node_id] < k){
    //                             k_list.push_back(edges[ee].node_id);
    //                             visited[edges[ee].node_id] = 1;
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //         for(auto& n : left[head]){
    //             if(!visited[n]){
    //                 deg[n]--;
    //                 if(deg[n] < k){
    //                     k_list.push_back(n);
    //                     visited[n] = 1;
    //                 }
    //             }
    //         }
    //     }
    //     // for(auto& p: redundant_edges){
    //     //    count++;
    //     // }
    // }
    // end = clock();
    // printf("Weighted set cover k-core time:%f s count : %d \n", ((double)(end - start) / (CLOCKS_PER_SEC )), count);

    // // // bfs
    // start = clock();
    // iterations =10;
    // count = 0;
    // while(iterations--){
    //     vector<int> visited(LM, 0);
    //     queue<int> q;
    //     for(int i = 0;i < LM;i++){
    //         if(!visited[i]){
    //             // 
    //             q.push(i);
    //             visited[i] = 1;
    //             while(!q.empty()){
    //                 int head = q.front();
    //                 q.pop();
    //                 // traverse the neighbors
    //                 for(int e = Vnode[head+input.size()];e != -1; e = next[e]){
    //                     if(IsEdge[e]){
    //                         for(int ee = Vnode[edges[e].node_id]; ee!= -1; ee = next[ee]){
    //                             if(IsEdge[ee] && edges[ee].isVirtual == 1 && !visited[edges[ee].node_id]){
    //                                 q.push(edges[ee].node_id);
    //                                 visited[edges[ee].node_id] = 1;
    //                             }
    //                         }
    //                     }        
    //                 }
    //                 for(auto& n : left[head]){
    //                     if(!visited[n]){
    //                         q.push(n);
    //                         visited[n] = 1;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    // end = clock();
    // printf("Weighted set cover bfs time:%f s count : %d \n", ((double)(end - start) / (CLOCKS_PER_SEC )), count);


    // // pagerank
    // // count the degree for each vertex.
    // start = clock();
    // // int iterations =10;
    // count = 0;
    // vector<int> deg(LM, 0);
    // for(int i = input.size(), s = input.size() + LM;i < s;i++){
    //     count = 0;
    //     for(int e = Vnode[i];e != -1; e = next[e]){
    //         if(IsEdge[e]){
    //             for(int ee = Vnode[edges[e].node_id]; ee!= -1; ee = next[ee]){
    //                 if(IsEdge[ee] && edges[ee].isVirtual == 1){
    //                     count++;
    //                 }
    //             }
    //         }
    //     }
    //     deg[i - s] = count;
    // }
    // for(int i = 0;i < LM;i++){
    //     deg[i] += left[i].size();
    // }

    // // perform iterations.
    // float* prev = new float[LM];
    // float* curr = new float[LM];
    // float damp = 0.85f;
    // // start = clock();
    // iterations = 1;
    // count = 0;
    // while(iterations--){
    //     for(int i = 0;i < LM;i++){
    //         curr[i] = (1.0f - damp) / LM;
    //         float t = 0.0f;
    //         for(int e = Vnode[i+input.size()];e != -1; e = next[e]){
    //             if(IsEdge[e]){
    //                 for(int ee = Vnode[edges[e].node_id]; ee!= -1; ee = next[ee]){
    //                     if(IsEdge[ee] && edges[ee].isVirtual == 1){
    //                         t += prev[edges[ee].node_id] / deg[edges[ee].node_id];
    //                     }
    //                 }
    //             }        
    //         }
    //         for(auto& n : left[i]){
    //             t += prev[n] / deg[n];
    //         }
    //         curr[i] += t * damp;
    //     }
    //     swap(prev, curr);
    // }
    
    // end = clock();
    // printf("Weighted set cover pagerank time:%f s count : %d \n", ((double)(end - start) / (CLOCKS_PER_SEC )), count);    
    // delete[] prev;
    // delete[] curr;

    delete[] cost;
    delete[] IsEdge;
    delete[] src;
    delete[] left_tag;
    delete[] right_tag;
    delete[] Vnode;
    delete[] counter;
}


// advanced search, combine set cover and search method.
void GraphDeduplicator::deduplicateWithNegaEdge(){
    printf("****************Negative Edge enhancement****************\n");

    ui v_num = input.size() + RM + LM;
    Vnode = new ui[v_num];
    memset(Vnode, -1, sizeof(ui) * v_num);
    // construct the threepart graph
    for(ui s = input.size(), i = 0;i < s;i++){
        NODE& vn = input[i];
        for(auto& ln : vn.left){
            insert_edge(i, ln, 1);
        }
        for(auto& rn : vn.right){
            insert_edge(i, rn, 0);
        }
    }
    
    set<pair<pair<int, int>, int>> redundant_edges;
    set<pair<pair<int, int>, int>> negative_edges;
    // record weight for each edge.
    ui* counter = new ui[edges.size()];
    memset(counter, 0, sizeof(ui) * edges.size());

    ui* left_tag = new ui[LM];
    ui* right_tag = new ui[RM];
    memset(left_tag, 0, sizeof(ui) * LM);
    memset(right_tag, 0, sizeof(ui) * RM);

    for(int i = 0, s = input.size();i < s;i++){
        NODE& node_i = input[i];
        for(int j = i + 1;j < s;j++){
            NODE& node_j = input[j];
            if (*(node_j.left.begin()) > *(--node_i.left.end())) {
                break;
            }
            
            // detection for conflict can be improved.
            set<int> node_i_j_left = findCommonNeighbour(node_i.left, node_j.left);
            if (node_i_j_left.size() == 0) continue;
            set<int> node_i_j_right = findCommonNeighbour(node_i.right, node_j.right);
            if (node_i_j_right.size() == 0) continue;

            // count
            for(auto& l : node_i_j_left){
                left_tag[l] = 1;
            }
            for(auto& r : node_i_j_right){
                right_tag[r] = 1;
            }

            // traverse the adjacent list of vertex i and j.
            for(int e = Vnode[i]; e != -1; e = next[e]){
                if((edges[e].isVirtual == -1 && left_tag[edges[e].node_id] == 1 )){
                    counter[e] += node_i_j_right.size() ;
                } 
                else if(edges[e].isVirtual == 1 &&right_tag[edges[e].node_id] == 1 ){
                    counter[e] += node_i_j_left.size() ;
                }
            }
            
            for(int e = Vnode[j]; e != -1; e = next[e]){
                if((edges[e].isVirtual == -1 && left_tag[edges[e].node_id] == 1 )){
                    counter[e] += node_i_j_right.size() ;
                } 
                else if(edges[e].isVirtual == 1 &&right_tag[edges[e].node_id] == 1 ){
                    counter[e] += node_i_j_left.size() ;
                }
            }

            // recover the state
            for(auto& l : node_i_j_left){
                left_tag[l] = 0;
            }
            for(auto& r : node_i_j_right){
                right_tag[r] = 0;
            }
        }
    }

    ui conflicts = 0;
    for(int i = 0;i <edges.size();i++){
        conflicts += counter[i];
    }
    conflicts /= 4;
    printf("there are %d conflicts.\n", conflicts);
    
    ui* src = new ui[edges.size()];
    for(int i = 0, s = input.size() + RM + LM;i < s;i++){
        for(int e = Vnode[i];e != -1;e = next[e]){
            src[e] = i;
        }
    }

    ui max_freq = 0;
    for(int i = 0;i < edges.size();i++){
        max_freq = max(max_freq, counter[i]);
    }
    // delete the edge with highest score until no conflicts exist.
    
    ui* Order = new ui[max_freq + 1];
    memset(Order, -1, sizeof(ui) * (max_freq+1));
    ui* Onext = new ui[edges.size()];
    bool* IsEdge = new bool[edges.size()];
    memset(IsEdge, -1, sizeof(bool) * edges.size());

    for(int i = 0;i < input.size();i++){
        for(int e = Vnode[i]; e != -1; e = next[e]){
            Onext[e] = Order[counter[e]];
            Order[counter[e]] = e;
        }
    }
    
    int delete_edges = 0;
    ui times = 0;
    ui* tag = nullptr;
    // resolve the conflicts until empty.
    for(int i = max_freq;i > 0;){
        while(i > 0 && Order[i] == -1){ i--;}
        if(i <= 0) break;
        int j = -1;
        // printf("i:%d \n", i);
        for(j = Order[i]; j != -1;){
            int tmp = Onext[j];
            if(IsEdge[j]) {
                if(counter[j] < i) {
                    Onext[j] = Order[counter[j]];
                    Order[counter[j]] = j;
                }
                else{
                    IsEdge[j] = 0;
                    // if(times <= 0) exit(1);
                    // printf("there are %d conflicts, %d left after this round.\n",conflicts, conflicts - counter[j]);
                    conflicts -= counter[j];
                    times--;
                    counter[j] = 0;
                    int s = src[j],  d = edges[j].node_id, type = edges[j].isVirtual;
                    
                    delete_edges++;
                    if(type == -1){    // d on the left side
                        // printf("delete (%d, %d)L\n", s, d);
                        d += input.size();
                        tag = right_tag;
                    }
                    else if(type == 1){     //d on the right side
                        // printf("delete (%d, %d)R\n", s, d);
                        d += input.size() + LM;
                        tag = left_tag;
                    }

                    for(int e = Vnode[d]; e != -1;e = next[e]){
                        if(edges[e].isVirtual == 0 && edges[e].node_id == s){
                            IsEdge[e] = 0;
                            break;
                        }
                    }

                    for(int e = Vnode[s]; e != -1 ;e = next[e]){
                        if(edges[e].isVirtual == -type && IsEdge[e]){
                            tag[edges[e].node_id] = 1;
                        }
                    }
                    for(int e = Vnode[d]; e != -1;e = next[e]){
                        // check if this edge is valid.
                        // if it is virtual node.
                        if(IsEdge[e] && edges[e].isVirtual == 0 && edges[e].node_id != s){
                            int c = 0, t = -1;
                            for(int ve = Vnode[edges[e].node_id]; ve != -1;ve = next[ve]){
                                if(edges[ve].isVirtual == -type && IsEdge[ve] && tag[edges[ve].node_id] >= 1){
                                    counter[ve]--;
                                    tag[edges[ve].node_id]++;
                                    c++;
                                }
                                else if (edges[ve].isVirtual == type){
                                    if(type == -1 && (d - input.size() == edges[ve].node_id)){
                                        t = ve;
                                    }
                                    else if(type == 1 && (d - input.size() - LM == edges[ve].node_id)){
                                        t = ve;
                                    }
                                }
                            }
                            counter[t] -= c;
                        }
                    }
                    int c1 = 0, c2 = 0;
                    for(int e = Vnode[s]; e != -1;e = next[e]){
                        if(edges[e].isVirtual == -type){
                            if(IsEdge[e]){
                                // add the broken connections.
                                if(tag[edges[e].node_id] == 1){
                                    c1++;
                                    // redundant_edges.insert(make_pair(make_pair(d, edges[e].node_id), type));
                                }
                                else if(tag[edges[e].node_id] > 1){
                                    c2++;
                                }
                                counter[e] = counter[e] - tag[edges[e].node_id] + 1;
                                // tag[edges[e].node_id] = 0;
                            }
                            else{
                                c2++;
                            }
                        }
                    }
                    for(int e = Vnode[s]; e != -1;e = next[e]){
                        if(edges[e].isVirtual == -type){
                            if(IsEdge[e]){
                                counter[e] = counter[e] - tag[edges[e].node_id] + 1;
                                if(tag[edges[e].node_id] == 1){
                                    // c1++;
                                    if(c2 + 1 >= c1){
                                        redundant_edges.insert(make_pair(make_pair(d, edges[e].node_id), type));
                                    }
                                }
                                else if(tag[edges[e].node_id] > 1 && c2 + 1< c1){
                                    negative_edges.insert(make_pair(make_pair(j, edges[e].node_id), type));
                                }
                                
                                tag[edges[e].node_id] = 0;
                            }
                            // add the broken connections.
                            else{
                                if(c2 + 1 < c1){
                                    negative_edges.insert(make_pair(make_pair(j, edges[e].node_id), type));
                                }
                            }
                        }
                    }
                    Order[i] = tmp;
                    break;
                }
            }
            j = tmp;
            // delete this edge and update the priority for each influenced edge.
        }
        if(j == -1) Order[i]= -1;
    }

    printf("delete edges:%d\n", delete_edges);
    // postprocess the degree one cases.
    for(int i = 0;i < input.size();i++){
        int left_nb = 0, right_nb = 0, l = 0, r = 0;
        for(int e = Vnode[i]; e != -1;e = next[e]){
            if(IsEdge[e]){
                if(edges[e].isVirtual == -1){
                    left_nb++;
                    l = edges[e].node_id;
                }
                else if(edges[e].isVirtual == 1){
                    right_nb++;
                    r = edges[e].node_id;
                }
            }
        }
        if(left_nb <= 1 ||  right_nb <= 1){
            // collapse the virtual i and add the direct edges.
            if(left_nb == 1){
                for(int e = Vnode[i]; e != -1;e = next[e]){
                    if(IsEdge[e] && edges[e].isVirtual == 1){
                        redundant_edges.insert(make_pair(make_pair(l + input.size(), edges[e].node_id), -1));
                    }
                }   
            }
            else if(right_nb == 1){
                for(int e = Vnode[i]; e != -1;e = next[e]){
                    if(IsEdge[e] && edges[e].isVirtual == -1){
                        redundant_edges.insert(make_pair(make_pair(r + input.size() + LM, edges[e].node_id), 1));
                    }
                }
            }
            for(int e = Vnode[i]; e != -1;e = next[e]){
                IsEdge[e] = 0;
            } 
        }
    }

    int c = 0;
    ui relations = 0;
    for(int i = 0;i < input.size();i++){
        int left_nb = 0, right_nb = 0;
        for(int e = Vnode[i]; e != -1;e = next[e]){
            if(IsEdge[e]){
                c++;
                if(edges[e].isVirtual == -1){
                    left_nb++;
                }
                else if(edges[e].isVirtual == 1){
                    right_nb++;
                }
            }
        }
        relations += left_nb * right_nb;
    }
    c += redundant_edges.size() + negative_edges.size();
    relations += redundant_edges.size();
    
    vector<vector<int> > Eadj(edges.size());
    for(auto& pp : negative_edges){
        int e = pp.first.first, v = pp.first.second, type = pp.second;
        // int s = src[e], d = edges[e].node_id; 
        Eadj[e].push_back(v);
    }

    int t = 0; 
    for(int i = 0, s = Eadj.size();i < s;i++){
        if(Eadj[i].size() == 0) continue;
        // t++;
        auto& l = Eadj[i];
        int u = src[i], v = edges[i].node_id, dire = edges[i].isVirtual; 
        assert(u < input.size());
        if(dire == -1){
            tag = right_tag;
        }
        else if(dire == 1){
            tag = left_tag;
        }
        for(auto& n : l){
            tag[n] = 1;
        }
        // 
        for(int e = Vnode[u]; e != -1;e = next[e]){
            if(edges[e].isVirtual == -dire && !tag[edges[e].node_id]){
                relations++;
            } 
        }

        for(auto& n : l){ 
            tag[n] = 0;
        }
    }
    
    c += t;

    printf("Number of edges after deduplication: %d\n", c);
    printf("Greedy Weighted Set cover Compression ratio: %f%\n", c * 100.0 / this->exp_size);
    printf("relations: %d\n", relations);

    // delete[] cost;
    delete[] IsEdge;
    delete[] src;
    delete[] left_tag;
    delete[] right_tag;
    delete[] Vnode;
    delete[] counter;
}

// advanced search, combine set cover and search method.
void GraphDeduplicator::WeightedDeduplicateWithNegaEdge(){
    printf("****************Weighted Negative Edge enhancement****************\n");

    ui v_num = input.size() + RM + LM;
    Vnode = new ui[v_num];
    memset(Vnode, -1, sizeof(ui) * v_num);
    // construct the threepart graph
    for(ui s = input.size(), i = 0;i < s;i++){
        NODE& vn = input[i];
        for(auto& ln : vn.left){
            insert_edge(i, ln, 1);
        }
        for(auto& rn : vn.right){
            insert_edge(i, rn, 0);
        }
    }
    
    set<pair<pair<int, int>, int>> redundant_edges;
    set<pair<pair<int, int>, int>> negative_edges;
    // record weight for each edge.
    ui* counter = new ui[edges.size()];
    memset(counter, 0, sizeof(ui) * edges.size());

    ui* cost = new ui[edges.size()];
    memset(cost, 0, sizeof(ui) * edges.size());

    ui* left_tag = new ui[LM];
    ui* right_tag = new ui[RM];
    memset(left_tag, 0, sizeof(ui) * LM);
    memset(right_tag, 0, sizeof(ui) * RM);

    for(int i = 0, s = input.size();i < s;i++){
        NODE& node_i = input[i];
        for(int j = i + 1;j < s;j++){
            NODE& node_j = input[j];
            if (*(node_j.left.begin()) > *(--node_i.left.end())) {
                break;
            }
            
            // detection for conflict can be improved.
            set<int> node_i_j_left = findCommonNeighbour(node_i.left, node_j.left);
            if (node_i_j_left.size() == 0) continue;
            set<int> node_i_j_right = findCommonNeighbour(node_i.right, node_j.right);
            if (node_i_j_right.size() == 0) continue;

            // count
            for(auto& l : node_i_j_left){
                left_tag[l] = 1;
            }
            for(auto& r : node_i_j_right){
                right_tag[r] = 1;
            }

            // traverse the adjacent list of vertex i and j.
            for(int e = Vnode[i]; e != -1; e = next[e]){
                if((edges[e].isVirtual == -1 && left_tag[edges[e].node_id] == 1 )){
                    counter[e] += node_i_j_right.size() ;
                } 
                else if(edges[e].isVirtual == 1 && right_tag[edges[e].node_id] == 1 ){
                    counter[e] += node_i_j_left.size() ;
                }
            }
            
            for(int e = Vnode[j]; e != -1; e = next[e]){
                if((edges[e].isVirtual == -1 && left_tag[edges[e].node_id] == 1 )){
                    counter[e] += node_i_j_right.size() ;
                } 
                else if(edges[e].isVirtual == 1 && right_tag[edges[e].node_id] == 1 ){
                    counter[e] += node_i_j_left.size() ;
                }
            }

            // recover the state
            for(auto& l : node_i_j_left){
                left_tag[l] = 0;
            }
            for(auto& r : node_i_j_right){
                right_tag[r] = 0;
            }
        }
    }

    #ifdef DEBUG
    ui conflicts = 0;
    for(int i = 0;i < edges.size();i++){
        conflicts += counter[i];
    }
    assert(conflicts % 4 == 0);
    conflicts /= 4;
    printf("there are %d conflicts.\n", conflicts);
    #endif

    ui* src = new ui[edges.size()];
    for(int i = 0, s = input.size() + RM + LM;i < s;i++){
        for(int e = Vnode[i];e != -1;e = next[e]){
            src[e] = i;
        }
    }

    for(int i = input.size(), k = input.size() + LM ;i < k;i++){
        for(int e = Vnode[i];e != -1;e = next[e]){
            ui mid = edges[e].node_id;
            // tag the vertices visited on right side.
            for(int ve = Vnode[mid]; ve != -1;ve = next[ve]){
                if(edges[ve].isVirtual == 1){
                    right_tag[edges[ve].node_id] += 1;
                }   
            }
        }
        // tag the vertices visited on right side.
        
        for(int e = Vnode[i];e != -1;e = next[e]){
            ui mid = edges[e].node_id;
            ui c = 0;
            int t = -1;
            for(int ve = Vnode[mid]; ve != -1;ve = next[ve]){
                if(edges[ve].isVirtual == 1 && right_tag[edges[ve].node_id] == 1){
                    if(counter[ve] > 0){
                        cost[ve]++;
                    }
                    c++;
                } 
                else if(edges[ve].isVirtual == -1 && edges[ve].node_id == i - input.size()){
                    t = ve;
                }
            }   
            if(counter[t] > 0){
                cost[t] = c;
            }
        }
        for(int e = Vnode[i];e != -1;e = next[e]){
            ui mid = edges[e].node_id;
            // tag the vertices visited on right side.
            for(int ve = Vnode[mid]; ve != -1;ve = next[ve]){
                if(edges[ve].isVirtual == 1){
                    right_tag[edges[ve].node_id] = 0;
                }   
            }
        }
    }

    priority_queue<pair<float, int>> pq;
    // ui* tag = nullptr;
    for(int i = 0;i < edges.size();i++){
        if(counter[i] > 0){
            // printf("%d %d %d, %d\n", src[i], edges[i].node_id, edges[i].isVirtual, cost[i]);
            pq.push(make_pair( cost[i] * -1.0 / counter[i] , i));
        }
    }

    ui max_freq = 0;
    for(int i = 0;i < edges.size();i++){
        max_freq = max(max_freq, counter[i]);
    }
    // delete the edge with highest score until no conflicts exist.


    bool* IsEdge = new bool[edges.size()];
    memset(IsEdge, -1, sizeof(bool) * edges.size());

    int delete_edges = 0;
    ui* tag = nullptr;

    int cc = 0;
    // resolve the conflicts until empty.
    while(!pq.empty()){
        if(conflicts == 0) {
            break;
        } 
        pair<float, int> p = pq.top();
        // printf("%f \n", p.first);
        pq.pop();
        int j = p.second;
        // pop the edge with the highest score and update others' score.
        if(IsEdge[j] != 0 && counter[j] > 0 && fabsf(p.first - cost[j] * -1.0f / counter[j]) < EPSINON) {
            // printf("%f %f\n", p.first,  (cost[j] * -1.0 / counter[j]));
            IsEdge[j] = 0;
            conflicts -= counter[j];
            counter[j] = 0;
            int s = src[j],  d = edges[j].node_id, type = edges[j].isVirtual;
            
            delete_edges++;
            if(type == -1){    // d on the left side
                // printf("delete (%d, %d)L\n", s, d);
                d += input.size();
                tag = right_tag;
            }
            else if(type == 1){     //d on the right side
                // printf("delete (%d, %d)R\n", s, d);
                d += input.size() + LM;
                tag = left_tag;
            }

            for(int e = Vnode[d]; e != -1;e = next[e]){
                if(edges[e].isVirtual == 0 && edges[e].node_id == s){
                    IsEdge[e] = 0;
                    break;
                }
            }

            for(int e = Vnode[s]; e != -1 ;e = next[e]){
                if(edges[e].isVirtual == -type && IsEdge[e]){
                    tag[edges[e].node_id] = 1;
                }
            }
            for(int e = Vnode[d]; e != -1;e = next[e]){
                // check if this edge is valid.
                // if it is virtual node.
                if(IsEdge[e] && edges[e].isVirtual == 0 && edges[e].node_id != s){
                    int c = 0, t = -1;
                    for(int ve = Vnode[edges[e].node_id]; ve != -1;ve = next[ve]){
                        if(edges[ve].isVirtual == -type && IsEdge[ve] && tag[edges[ve].node_id] >= 1){
                            counter[ve]--;
                            tag[edges[ve].node_id]++;
                            c++;
                        }
                        else if (edges[ve].isVirtual == type){
                            if(type == -1 && (d - input.size() == edges[ve].node_id)){
                                t = ve;
                            }
                            else if(type == 1 && (d - input.size() - LM == edges[ve].node_id)){
                                t = ve;
                            }
                        }
                    }
                    counter[t] -= c;
                }
            }
            for(int e = Vnode[d]; e != -1;e = next[e]){
                // if it is virtual node.
                if(IsEdge[e] && edges[e].isVirtual == 0 && edges[e].node_id != s){
                    int c = 0, t = -1;
                    for(int ve = Vnode[edges[e].node_id]; ve != -1;ve = next[ve]){
                        if(edges[ve].isVirtual == -type && IsEdge[ve] && tag[edges[ve].node_id] >= 1){
                            if(tag[edges[ve].node_id] == 2) {
                                cost[ve]++;  
                                c++; 
                            }
                            if(counter[ve] > 0){
                                pq.push(make_pair(cost[ve] * -1.0f / counter[ve] , ve));
                            }
                        }
                        else if (edges[ve].isVirtual == type){
                            if(type == -1 && (d - input.size() == edges[ve].node_id)){
                                t = ve;
                            }
                            else if(type == 1 && (d - input.size() - LM == edges[ve].node_id)){
                                t = ve;
                            }
                        }
                    }
                    cost[t] += c;
                    if(counter[t] > 0){
                        pq.push(make_pair(cost[t] * -1.0f / counter[t], t));
                    }
                }
            }

            int c1 = 0, c2= 0;
            for(int e = Vnode[s]; e != -1;e = next[e]){
                if(edges[e].isVirtual == -type){
                    // add the broken connections.
                    // counter[e] = counter[e] - tag[edges[e].node_id] + 1;
                    if(IsEdge[e]){
                        if(tag[edges[e].node_id] == 1){
                        // cost[e]--;
                            c1++;
                            // redundant_edges.insert(make_pair(make_pair(d, edges[e].node_id), type));
                        }
                        else if(tag[edges[e].node_id] > 1){
                            c2++;
                        }
                    }
                    else{
                        c2++;
                    }
                    // if(counter[e] > 0){
                    //     pq.push(make_pair(cost[e] * -1.0f / counter[e], e));
                    // }
                    // tag[edges[e].node_id] = 0;
                }
            }
            
            for(int e = Vnode[s]; e != -1;e = next[e]){
                if(edges[e].isVirtual == -type){
                    if(IsEdge[e]){
                        counter[e] = counter[e] - tag[edges[e].node_id] + 1;
                        if(tag[edges[e].node_id] == 1){
                            cost[e]--;
                            // c1++;
                            if(c2 + 1 >= c1){
                                redundant_edges.insert(make_pair(make_pair(d, edges[e].node_id), type));
                            }
                        }
                        else if(tag[edges[e].node_id] > 1 && c2 + 1< c1){
                            negative_edges.insert(make_pair(make_pair(j, edges[e].node_id), type));
                        }
                        
                        if(counter[e] > 0){
                            pq.push(make_pair(cost[e] * -1.0f / counter[e], e));
                        }
                        tag[edges[e].node_id] = 0;
                    }
                    // add the broken connections.
                    else{
                        if(c2 + 1 < c1){
                            negative_edges.insert(make_pair(make_pair(j, edges[e].node_id), type));
                        }
                    }
                }
                
            }
        }
    }

    // printf("cc: %d\n", cc);

    for(int i = 0;i < edges.size();i++){
        if(counter[i] > 0){
            printf("error\n");
        }
    }

    printf("delete edges:%d\n", delete_edges);

    // postprocess the degree one cases.
    for(int i = 0;i < input.size();i++){
        int left_nb = 0, right_nb = 0, l = 0, r = 0;
        for(int e = Vnode[i]; e != -1;e = next[e]){
            if(IsEdge[e]){
                if(edges[e].isVirtual == -1){
                    left_nb++;
                    l = edges[e].node_id;
                }
                else if(edges[e].isVirtual == 1){
                    right_nb++;
                    r = edges[e].node_id;
                }
            }
        }
        if(left_nb <= 1 ||  right_nb <= 1){
            // collapse the virtual i and add the direct edges.
            if(left_nb == 1){
                for(int e = Vnode[i]; e != -1;e = next[e]){
                    if(IsEdge[e] && edges[e].isVirtual == 1){
                        redundant_edges.insert(make_pair(make_pair(l + input.size(), edges[e].node_id), -1));
                    }
                }   
            }
            else if(right_nb == 1){
                for(int e = Vnode[i]; e != -1;e = next[e]){
                    if(IsEdge[e] && edges[e].isVirtual == -1){
                        redundant_edges.insert(make_pair(make_pair(r + input.size() + LM, edges[e].node_id), 1));
                    }
                }
            }
            for(int e = Vnode[i]; e != -1;e = next[e]){
                IsEdge[e] = 0;
            } 
        }
    }

    int c = 0;
    ui relations = 0;
    for(int i = 0;i < input.size();i++){
        int left_nb = 0, right_nb = 0;
        for(int e = Vnode[i]; e != -1;e = next[e]){
            if(IsEdge[e]){
                c++;
                if(edges[e].isVirtual == -1){
                    left_nb++;
                }
                else if(edges[e].isVirtual == 1){
                    right_nb++;
                }
            }
        }
        relations += left_nb * right_nb;
    }
    c += redundant_edges.size() + negative_edges.size();
    relations += redundant_edges.size();
    
    vector<vector<int> > Eadj(edges.size());
    for(auto& pp : negative_edges){
        int e = pp.first.first, v = pp.first.second, type = pp.second;
        // int s = src[e], d = edges[e].node_id; 
        Eadj[e].push_back(v);
    }

    int t = 0; 
    for(int i = 0, s = Eadj.size();i < s;i++){
        if(Eadj[i].size() == 0) continue;
        // t++;
        auto& l = Eadj[i];
        int u = src[i], v = edges[i].node_id, dire = edges[i].isVirtual; 
        assert(u < input.size());
        if(dire == -1){
            tag = right_tag;
        }
        else if(dire == 1){
            tag = left_tag;
        }
        for(auto& n : l){
            tag[n] = 1;
        }
        // 
        for(int e = Vnode[u]; e != -1;e = next[e]){
            if(edges[e].isVirtual == -dire && !tag[edges[e].node_id]){
                relations++;
            } 
        }

        for(auto& n : l){ 
            tag[n] = 0;
        }
    }
    
    c += t;

    printf("Number of edges after deduplication: %d\n", c);
    printf("Greedy Weighted Set cover Compression ratio: %f%\n", c * 100.0 / this->exp_size);
    printf("relations: %d\n", relations);

    delete[] cost;
    delete[] IsEdge;
    delete[] src;
    delete[] left_tag;
    delete[] right_tag;
    delete[] Vnode;
    delete[] counter;
}

// count the number of edges after deduplication and report the compression ratio.
void GraphDeduplicator::report_result(){
    unsigned int counter = 0;
    for(int s = out_left.size(), i = 0;i < s;i++){
        vector<OUT_NODE>& x = out_left[i];
        for(auto& n : x){
            if(n.isVirtual){
                counter++;
            }  
        }
    }

    for(int s = out_right.size(), i = 0;i < s;i++){
        counter += out_right[i].size();
    }

    printf("Number of edges after deduplication: %d \n", counter);
    printf("Compression ratio: %f%\n", 100.0f * counter / this->exp_size);


    ui relations = 0;
    for(auto& mid : input){
        relations += ((mid.left.size() - mid.left_del.size()) * (mid.right.size() - mid.right_del.size()));
    }

    for(int s = out_left.size(), i = 0;i < s;i++){
        vector<OUT_NODE>& x = out_left[i];
        for(auto& n : x){
            if(!n.isVirtual){
                relations++;
            }  
        }
    }
    printf("relations :%d\n", relations);
    // assert(relations == this->exp_size);
}


void GraphDeduplicator::print_graph(){
    printf("***************Left side*******************\n");
    for(int s = out_left.size(), i = 0;i < s;i++){
        vector<OUT_NODE>& x = out_left[i];
        for(auto& n : x){
            if(n.isVirtual){
                printf("%d %d MID \n", i, input[n.node_id].input_node_id);
            }
            else{
                printf("%d %d R\n", i, n.node_id);
            }
        }
    }

    printf("***************Right side*******************\n");
    for(int s = out_right.size(), i = 0;i < s;i++){
        vector<OUT_NODE>& x = out_right[i];
        for(auto& n : x){
            if(n.isVirtual){
                printf("%d %d MID \n", i, input[n.node_id].input_node_id);
            }
            else{
                printf("%d %d L\n", i, n.node_id);
            }
        }
    }
}

vector<vector<int>> GraphDeduplicator::get_neighbour() {
    return neighbour;
}

vector<int> GraphDeduplicator::get_weight() {
    return weight;
}

// pair<VVON, VVON> GraphDeduplicator::get_result() {
//     return make_pair(out_left, out_right);
// }

