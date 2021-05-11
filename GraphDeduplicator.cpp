#include <cstring>
#include <fstream>
#include <set>
#include <map>
#include <vector>
#include <queue>
#include <algorithm>
#include <iterator>
#include <utility>
#include <iostream>
#include <ctime>
#include "GraphDeduplicator.h"
#include "assert.h"


using namespace std;

vector<string> split_line(string line, string delimiter = "\r\t\n ") {
    std::size_t pre_pos, pos;
    vector<string> res;
    pre_pos = 0;
    pos = line.find_first_of(delimiter, pre_pos);
    while (pos != std::string::npos) {
        res.push_back(line.substr(pre_pos, pos - pre_pos));
        pre_pos = pos + 1;
        pos = line.find_first_of(delimiter, pre_pos);
    }
    if (pre_pos != pos) {
        res.push_back(line.substr(pre_pos, pos - pre_pos));
    }
    return res;
}


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

GraphDeduplicator::GraphDeduplicator(char* filename) {
    // read file
    ifstream ifs(filename, std::ifstream::in);
    if (ifs) {
        string line;
        while (getline(ifs, line)) {
            ++N;
            if (line[0] == '#' || line[0] == '%') {
                continue;
            }
            vector<string> content = split_line(line);
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
            raw_edge_tot += (left_len + right_len);
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

bool compare(const NODE& x, const NODE& y) {
    int x_min = *(x.left.begin());
    int y_min = *(y.left.begin());
    if (x_min < y_min) {
        return true;
    } else if (x_min == y_min) {
        return *(--x.left.end()) < *(--y.left.end()); //  
    }
    return false;
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
    // sort
    sort(input.begin(), input.end(), compare); 
    
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
                    counter[j] = 0;
                    int s = src[j],  d = edges[j].node_id, type = edges[j].isVirtual;
                    ui* tag = nullptr;
                    if(type == -1){    // d on the left side
                        d += input.size();
                        tag = right_tag;
                    }
                    else if(type == 1){     //d on the right side
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
    // print the direct edges.
    // for(auto& p : redundant_edges){
    //     auto& pp = p.first;
    //     int type = p.second;
    //     if(type == -1){
    //         printf("L %d -> %d\n",pp.first - input.size(), pp.second);
    //     }
    //     else if(type == 1){
    //         printf("R %d -> %d\n",pp.first - input.size() - LM, pp.second);
    //     }
    // }

    printf("****************Set Cover****************\n");
    printf("Number of edges after deduplication: %d\n", c);
    printf("Set cover Compression ratio: %f%\n", c * 100.0 / this->exp_size);
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
    sort(input.begin(), input.end(), compare); 

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

            // int num_exp_ij = node_i.left.size() * node_i.right.size() + node_j.left.size() * node_j.right.size() 
            //                                          - node_i_j_left.size() * node_i_j_right.size();
            
            // int c = node_i.left.size() + node_j.left.size() + node_j.right.size() + node_j.left.size();


            // weight[i] += max(max((int)(num_exp_ij - c + node_i_j_left.size() - node_i_j_left.size() * (node_j.right.size() - node_i_j_right.size())), 1), 
            //                     max((int)(num_exp_ij - c + node_i_j_right.size() - node_i_j_right.size() * (node_j.left.size() - node_i_j_left.size())), 1));

            // weight[j] += max(max((int)(num_exp_ij - c + node_i_j_left.size() - node_i_j_left.size() * (node_i.right.size() - node_i_j_right.size())), 1), 
            //                     max((int)(num_exp_ij - c + node_i_j_right.size() - node_i_j_right.size() * (node_i.left.size() - node_i_j_left.size())), 1));
        }
    }
}

void GraphDeduplicator::compute_mis_benefit(vector<int>& mis){
    int counter = 0;
    for(auto& v : mis){
        counter += input[v].left.size() * input[v].right.size();
    }
    printf("Edges covered by this mis: %d\n", counter);
}

void GraphDeduplicator::assign_vertex_weight(int k) {
    clock_t start = clock();
    weight.resize(N, 1);
    // if (k == 1) {
    //     weight.clear();
    //     weight.resize(N, 1);
    // }
    // else if (k == 2) {
    //     for (int i = 0; i < N; ++i) {
    //         int left = input[i].left.size();
    //         int right = input[i].right.size();
    //         // int id = input[i].node_id;
    //         // weight < 0 ?
    //         weight[i] = max(left * right - (left + right), 1);
    //         // weight[i] = left * right;
    //     }
    // }
    // else if (k == 3) {
        // for (int i = 0; i < N; ++i) {
        //     for(auto& n : neighbour[i]){

        //     }
        //     // int id = input[i].node_id;
            
        //     // --weight[id];
        //     // vector<NODE_COMMON_NEIGHBOUR> cn = input[i].common_neighbour;
        //     // for (int j = 0; j < cn.size(); ++j) {
        //     //     weight[id] += (cn[j].left_common_neighbour.size() + cn[j].right_common_neighbour.size()); 
        //     // }
        // }
    // }
    // clock_t end = clock();
    // printf("Assign weight time:%fs\n", (double)(end-start) / CLOCKS_PER_SEC);
}

void GraphDeduplicator::build_conflict_graph() {
    construct_graphTopology();
    assign_vertex_weight(1);
}

// greedy real node first deduplication
void GraphDeduplicator::dedup1(){
    build_revertedIDX_for_realNodes();
    set<pair<int,int>> redundant_edges;
    for(auto& real_node : out_left){
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
    
    // clean up 
    for(auto& n : input){
        n.left_del.clear();
        n.right_del.clear();
    }
    out_left.clear();
    out_right.clear();
}

int GraphDeduplicator::count_expand_edges(){
    vector<vector<int>> tmp_left;
    tmp_left.resize(LM + 1, {});
    for(int s = input.size(), i = 0;i < s;i++){
        NODE& node = input[i];
        for(auto& Lrn : node.left){
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
    for (int i = 0; i <= LM; ++i) {
        out_left.push_back({});
    }
    for (int i = 0; i <= RM; ++i) {
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
    vector<int> unconflict_set;
    unconflict_set.resize(input.size(), 0);
    set<pair<int,int>> redundant_edges;

    // label each vertex which is mutually unconflict.
    for(auto n : mis){
        unconflict_set[n] = 1;
    }

    // deduplicate one by one
    for(auto id : mvc){
        for(auto& ngb : neighbour[id]){
            // conflict exist, resolve it
            if(unconflict_set[ngb] == 1){
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

    delete[] tag;
    printf("****************Search****************\n");
    report_result();
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

pair<VVON, VVON> GraphDeduplicator::get_result() {
    return make_pair(out_left, out_right);
}