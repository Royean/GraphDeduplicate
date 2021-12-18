#pragma once
#include "Bi2Tri.h"

bool degree_compar(const pair<int,int>& a, const pair<int,int>& b){
    if(a.first < b.first) return true;
    return false;
}

bool subgraph_size_comp(){

}

// power-law distribution.
Bi2Tri::Bi2Tri(float alpha, int n, int m){
    left.resize(n);
    right.resize(m);

    vector<float> deg;
    deg.resize(m,1.0);

    float sum = m;

    vector<int> L;
    for(int i = 0;i < m;i++){
        L.push_back(i);
        // left[0].push_back(i);
    }

    // int e = 0;
    srand(unsigned(time(0)));
    for(int i = 0;i < n;i++){
        for(int j = 0;j < m;j++){
            float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            if(r < deg[j] * 1.0 / sum){
                left[i].push_back(L[j]);
                right[L[j]].push_back(i);
                e++;
                deg[j] += 1.0 * alpha;
                sum += 1.0 * alpha;
            }
        } 
        // random_shuffle(L.begin(), L.end(), Rand);
        // for(int j = 0;j < (int)( * m);j++){
            // left[i].push_back(L[j]);
            // right[L[j]].push_back(i);
            // e++;
        // }
    }
    printf("bigraph reading done......\n");
    printf("left nodes: %ld, right nodes: %ld, edges: %d\n", left.size(), right.size(), e);
}

// Bi2Tri::Bi2Tri(int n, int m, float prob){
//     // int n = n, m;
//     left.resize(n);
//     right.resize(m);
//     L_Deg.resize(n);
//     R_Deg.resize(m);
//     L_valid.resize(n);
//     R_valid.resize(m);

//     vector<int> L;
    
//     for(int i = 0;i < m;i++){
//         L.push_back(i);
//     }

//     // int e = 0;
//     srand(unsigned(time(0)));
//     for(int i = 0;i < n;i++){
//         random_shuffle(L.begin(), L.end(), Rand);
//         for(int j = 0;j < (int)(prob * m);j++){
//             left[i].push_back(L[j]);
//             right[L[j]].push_back(i);
//             e++;
//         }
//     }
//     printf("bigraph reading done......\n");
//     printf("left nodes: %ld, right nodes: %ld, edges: %d\n", left.size(), right.size(), e);
// }


Bi2Tri::Bi2Tri(string filename){
    ifstream ifs("./bigraph/"+ filename, std::ifstream::in);
    int e = 0;
    if (ifs) {
        string line;
        getline(ifs, line);
        int n = atoi(line.c_str());
        left.resize(n);
        getline(ifs, line);
        int m = atoi(line.c_str());
        right.resize(m);
        this->ln = n;
        this->rn = m;
        L_Deg = new int[n];
        R_Deg = new int[m];
        L_valid = new int[n];
        R_valid = new int[m];
        while (getline(ifs, line)) {
            if (line[0] == '#' || line[0] == '%') {
                continue;
            }
            vector<string> content = split_line(line, "\r\t\n ");
            int l = atoi(content[0].c_str());
            int k = 1;
            int left_len = atoi(content[k++].c_str());
            // set<int> left_tmp;
            for (int i = 0; i < left_len; i++) {
                assert(k < content.size());
                int temp = atoi(content[k++].c_str());
                assert(l < left.size());
                left[l].push_back(temp);
                assert(temp < right.size());
                right[temp].push_back(l);
                e++;
            }
        }
    }
    ifs.close();
    this->e = e;

    for(int i = 0, s = left.size();i < s;i++){
        L_Deg[i] = left[i].size();
        max_degree = max(max_degree, L_Deg[i]);
        L_valid[i] = 1;
    }

    for(int i = 0, s = right.size();i < s;i++){
        R_Deg[i] = right[i].size();
        max_degree = max(max_degree, R_Deg[i]);
        R_valid[i] = 1;
    }

    printf("reading done......\n");
    printf("|L|: %ld, |R|: %ld, |E|: %d\n", left.size(), right.size(), e);
}

void Bi2Tri::subgraph_reduce(){

}

void Bi2Tri::degree_one_reduce(){
    vector<int> L_degree_one, R_degree_one;
    for(int i = 0, s = left.size();i < s;i++){
        if(left[i].size() == 1){
            L_degree_one.push_back(i);
        }
    }
    
    for(int i = 0, s = right.size();i < s;i++){
        if(right[i].size() == 1){
            R_degree_one.push_back(i);
        }
    }
    while(!L_degree_one.empty() || !R_degree_one.empty()){
        while (!L_degree_one.empty()){
            int u = L_degree_one.back();
            L_degree_one.pop_back();
            // assert(u < left.size());
            if(L_Deg[u] == 1){
                L_valid[u] = 0;
                for(auto r : left[u]) if(R_valid[r]){
                    // assert(r < right.size());
                    R_Deg[r]--;
                    if(R_Deg[r] == 1){
                        R_degree_one.push_back(r);
                    }
                }
            }
            
        }
        while(!R_degree_one.empty()){
            int u = R_degree_one.back();
            R_degree_one.pop_back();
            // if(u >= right.size()){
            //     // cerr << u << " " << right.size() << endl;
            //     // printf("%d %d", u, right.size());
            // }
            // assert(u < right.size());
            if(R_Deg[u] == 1){
                R_valid[u] = 0;
                for(auto l : right[u]) if (L_valid[l]){
                    // assert(l < left.size());
                    L_Deg[l]--;
                    if(L_Deg[l] == 1){
                        L_degree_one.push_back(l);
                    }
                }
            }
            
        }
    }
    printf("Degree one reduce done...\n");
}

void Bi2Tri::upperbound_reduce(){
    // vector<int> L_counter, R_counter;
    ui* L_counter = new ui[left.size()];
    ui* R_counter = new ui[right.size()];
    char* L_skip = new char[left.size()];
    char* R_skip = new char[right.size()];
    memset(L_counter, 0, sizeof(ui) * left.size());
    memset(R_counter, 0, sizeof(ui) * right.size());

    memset(L_skip, 0, sizeof(char) * left.size());
    memset(R_skip, 0, sizeof(char) * right.size());

    for(int l = 0, s = left.size();l < s; l++){
        if(L_valid[l] == 1 && L_skip[l] == 0){
            char ok = 1;
            vector<int> N;
            for(int i = 0, rs = left[l].size();i < rs;i++){
                int r = left[l][i];
                if(R_valid[r] == 0 || R_skip[r] == 1) continue;
                int two_counter = 0;
                for(int j = 0, lls = right[r].size();j < lls;j++){
                    int ll = right[r][j];
                    if(ll != l && L_valid[ll] == 1){
                        L_counter[ll]++;
                        if(L_counter[ll] == 1){
                            N.push_back(ll);
                        }
                        if(L_counter[ll] == 2){
                            two_counter++;
                            if(two_counter > 1){
                                ok = 0;
                                L_skip[ll] = 1;
                                R_skip[r] = 1;
                            }
                        }
                        if(L_counter[ll] == 3){
                            ok = 0;
                            // tag the one-hop vertices.
                            L_skip[ll] = 1;
                        }
                    }
                }
            }
            if(ok == 1){
                L_valid[l] = 0;
            }
            for(auto& x : N){
                L_counter[x] = 0;
            }
        }
    }

    for(int l = 0, s = right.size();l < s; l++){
        if(R_valid[l] == 1 && R_skip[l] == 0){
            char ok = 1;
            vector<int> N;
            for(int r = 0, rs = right[l].size();r < rs;r++){
                int k = right[l][r];
                if(L_valid[r] == 0) continue;
                int two_counter = 0;
                for(int ll = 0, lls = left[k].size();ll < lls;ll++){
                    int llk = left[k][ll];
                    if(llk != l && R_valid[llk] == 1){
                        R_counter[llk]++;
                        if(R_counter[llk] == 1){
                            N.push_back(llk);
                        }
                        if(R_counter[llk] == 2){
                            two_counter++;
                            if(two_counter > 1){
                                ok = 0;
                                R_skip[llk] = 1;
                            }
                        }
                        if(R_counter[llk] == 3){
                            ok = 0;
                            // tag the one-hop vertices.
                            R_skip[llk] = 1;
                        }
                    }
                }
            }
            if(ok == 1){
                R_valid[l] = 0;
            }
            for(auto& x : N){
                R_counter[x] = 0;
            }
        }
    }
    delete[] L_skip;
    delete[] R_skip;
    delete[] L_counter;
    delete[] R_counter;
    printf("two-hop reduce done...\n");
}

void Bi2Tri::reduce(){
    degree_one_reduce();
    summary();
    upperbound_reduce();
    summary();
    printf("Bigraph reduce done...\n");
}

vector<NODE> Bi2Tri::duplicate_allowed_compress(){
    int* tagL = new int[left.size()];
    memset(tagL, 0, sizeof(int) * left.size());

    int* tagR = new int[right.size()];
    memset(tagR, 0, sizeof(int) * right.size());

    ui*  ver_intersec = new ui[left.size()];
    memset(ver_intersec, 0, sizeof(ui) * left.size());
    vector<vector<int> > isEdge(left);
    
    int mid_count = 0;
    vector<NODE> input;
    int c = 0;
      
    vector<pair<int,int> > L;
    for(int i = 0;i < left.size();i++){
        if(L_valid[i] == 1){
            L.push_back(make_pair(left[i].size(), i));
        }
    }

    /**
     * @brief try different vertex ordering.
     *  1. degree-degeneracy order
     *  2. subgraph-degeneracy order
     */
    sort(L.begin(), L.end(), degree_compar);

    int* visited = new int[left.size()];
    memset(visited, 0, sizeof(ui) * left.size());

    int* jump = new int[left.size()];
    memset(jump, 0, sizeof(ui) * left.size());

    int* start = new int[left.size()];
    
    int* local_degree = new int[left.size()];

    int* local_dest = new int[right.size()];

    for(int q = 0, s = L.size();q < s;q++){ 
        int i = L[q].second;
        // if(jump[i] == 1){
        //     continue;
        // } 
        ui degs = 0;
        vector<int> one_hop, two_hop;
        tagL[i] = 1;
        one_hop.resize(max_degree);
        two_hop.resize(max_degree * max_degree);
        for(int j = 0;j < left[i].size();j++) {
            int r = left[i][j];
            tagR[r] = 1;
            one_hop.push_back(r);
            for(int k = 0;k < right[r].size();k++) {
                int l = right[r][k];
                if(l != i){
                    tagL[l]++;
                    if(tagL[l] == 2){      // // two_hop store the vertices with degree >= 2;
                        two_hop.push_back(l);
                    }
                }
            }
        }

        
        // apply the reductions for the extracted subgraph.
        subgraph_reduce();
        
        ui vs = left[i].size();
        priority_queue<pair<int, int>> pq;
        float threshold = 0;
        int delta_max = -1; // neg-infinity       
        // for each candidate, calculate the profit.
        for(auto& x : two_hop){
            int d_count = 0, v_count = 0;
            for(int j = 0;j < left[x].size();j++) {
            // for(auto& xr : left[x]) {
                int xr = left[x][j];
                if(tagR[xr] == 1){
                    v_count++;
                }
            }
            ver_intersec[x] = v_count;
            delta_max = max(delta_max, v_count);
            pq.push(make_pair(ver_intersec[x], x));   
        }

        // choose the vertex with highest score.
        vector<int> Ls;
        Ls.push_back(i);
        ui Rs = 0;
        int times = 1;
        while(!pq.empty()){
            // printf("%d\n", pq.size());
             // update the thresho ld.
            pair<int, int>  p = pq.top();
            pq.pop();
            if(p.first <= (int)floor(threshold)){
                break;
            }
            if(p.first != ver_intersec[p.second]){
                continue;
            }
            Ls.push_back(p.second);
            tagL[p.second] = 0;
            int temp = 0;
            for(int j = 0;j < left[p.second].size();j++){
            // for(auto& x : left[p.second]){ 
                int x = left[p.second][j];
                if(tagR[x] == times){
                    tagR[x]++;   
                    temp++;
                }
            }
            times++;
            Rs = temp;
            
            vector<int> update; 
            for(int j = 0;j < left[i].size();j++) {
            // for(auto& x : left[i]){
                int x = left[i][j];
                if(tagR[x] == times - 1){
                    // update the v_diff and d_diff
                    for(auto& xl : right[x]) if(xl != i){
                        ver_intersec[xl]--;
                        // deg_intersec[xl] -= right[x].size();
                        if(visited[xl] == 0){
                            update.push_back(xl);
                            visited[xl] = 1;
                        }
                        // float s = ver_intersec[xl];
                        // pq.push(make_pair(s, xl));
                        // v_pq.push(make_pair(ver_intersec[xl], xl));
                    }
                }
            }

            for(auto& xl : update){
                int s = ver_intersec[xl];
                pq.push(make_pair(s, xl));
                // v_pq.push(make_pair(ver_intersec[xl], xl));
                visited[xl] = 0;
            }

            // add a new virtual node controlled by the parameter *redundancy*.
            vector<int> Rsol;
            for(int j = 0;j < left[i].size();j++){
                int r = left[i][j];
                if(tagR[r] == times){
                    Rsol.push_back(r);  
                }
            }
            
            if(Rsol.size() > 1 && Ls.size() > 1){
                set<int> left_(Ls.begin(), Ls.end());
                set<int> right_(Rsol.begin(), Rsol.end());
                // printf("%d %d\n",right_.size() , left_.size());
                input.push_back(NODE(mid_count++, left_, right_));
            }
            threshold = Ls.size() * Rs * 1.0 / (Ls.size() + 1);
        }

        
        for(auto& x : two_hop){
            tagL[x] = 0;
        }
        tagL[i] = 0;
        vector<int> Rsol, Ridx;
        
        // pick out common neighbors of vertices in Ls.
        for(int j = 0;j < left[i].size();j++){
            int r = left[i][j];
            if(tagR[r] == times){
                Rsol.push_back(r);  
            }
        }
        
        // generate the virtual node.
        if(Rsol.size() > 1 && Ls.size() > 1){
            for(auto& ls : Ls){
                jump[ls] = 1;
            }
            set<int> left_(Ls.begin(), Ls.end());
            set<int> right_(Rsol.begin(), Rsol.end());
            // printf("%d %d\n",right_.size() , left_.size());
            input.push_back(NODE(mid_count++, left_, right_));
        }

        for(auto& r : left[i]){
            tagR[r] = 0;
        }
    }

    delete[] tagL;
    delete[] tagR;
    delete[] ver_intersec;
    delete[] visited;
    delete[] start;
    delete[] local_degree;
    delete[] local_dest;

    printf("%d mid nodes are generated...\n", mid_count);

    set<pair<int,int>> T;
    int ec = 0;
    int k = 0;
    for(auto& node : input){
        k += node.left.size() * node.right.size();
        ec += node.left.size() + node.right.size();
        for(auto& l : node.left){
            for(auto& r : node.right){
                T.insert(make_pair(l,r));
            }
        }
    }    
    printf("Trigrahph size:%d Expanded: %d direct edges:%d\n", ec, T.size(), e - T.size());
    // assert(e -T.size() == c);
    printf("comrepssion ratio: %f %\n", (e - T.size() + ec) * 100.0 / this->e);   
    return input;
}

// convert the bigraph to trigraph
// vector<NODE> Bi2Tri::convert(){
//     int* tagL = new int[left.size()];
//     memset(tagL, 0, sizeof(int) * left.size());

//     int* tagR = new int[right.size()];
//     memset(tagR, 0, sizeof(int) * right.size());

//     ui* deg_intersec = new ui[left.size()];
//     memset(deg_intersec, 0, sizeof(ui) * left.size());

//     ui* ver_intersec = new ui[left.size()];
//     memset(ver_intersec, 0, sizeof(ui) * left.size());

//     // ui* deg= new ui[right.size()];
//     // memset(deg, 0, sizeof(ui) * right.size());

//     int mid_count = 0;
//     // int midsize = 0;
//     vector<NODE> input;
    
//     int c = 0;
//     for(int i = 0, s = left.size();i < s;i++) {
//         // extract the two-hop neighbors for each vertex.
//         // printf("%d th vertex is processed...\n", i);
//         ui degs = 0;
//         vector<int> one_hop, two_hop;
//         tagL[i] = 1;
//         for(auto& r : left[i]){
//             degs += right[r].size();
//             tagR[r] = 1;
//             for(auto& l : right[r]){
//                 if(!tagL[l]){
//                     two_hop.push_back(l);
//                     tagL[l] = 1;
//                 }
//             }
//         }

//         ui vs = left[i].size();
//         priority_queue<pair<float, int>> pq;
//         priority_queue<pair<int, int>> v_pq;
//         float threshold = 0;

//         int delta_max = -1; // neg-infinity
//         // calculate the difference for each vertex
//         for(auto& x : two_hop){
//             int d_count = 0, v_count = 0;
//             for(auto& xr : left[x]){
//                 if(tagR[xr] == 1){
//                     v_count++;
//                     d_count += right[xr].size();
//                 }
//             }
//             deg_intersec[x] = d_count;
//             ver_intersec[x] = v_count;
//             delta_max = max(delta_max, v_count);
//             pq.push(make_pair(ver_intersec[x], x));   
//             // pq.push(make_pair(deg_intersec[x] * 1.0 / (ver_intersec[x] + 1), x));   
//             v_pq.push(make_pair(v_count,x));
//         }

//         // printf("%f\n", threshold);
//         vector<int> Ls;
//         Ls.push_back(i);
//         ui Rs = 0;
//         int times = 1;
//         while(!pq.empty()){
//             // printf("%d %d\n", delta_max, (int)floor(threshold));
//              // update the thresho ld.
//             if(delta_max <= (int)floor(threshold)){
//                 break;
//             }
//             pair<float, int>  p = pq.top();
//             pq.pop();
//             // assert(p.second < left.size());
//             if(fabsf(p.first - ver_intersec[p.second]) > EPSINON){
//                 continue;
//             }
//             // if(fabsf(p.first - deg_intersec[p.second] * 1.0 / (ver_intersec[p.second] + 1)) > EPSINON){
//             //     continue;
//             // }
//             // printf("merge %d\n", p.second);
//             Ls.push_back(p.second);
//             tagL[p.second] = 0;
//             int temp = 0;
//             for(auto& x : left[p.second]){ 
//                 // assert(x < right.size());
//                 if(tagR[x] == times){
//                     // printf("%d %d\n", x, tagR[x]);
//                     tagR[x]++;   
//                     temp++;
//                     // printf("%d %d\n", x, tagR[x]);
//                 }
//             }
//             // printf("\n");
//             // for(auto& r : left[i]){
//             //     printf("tagR:%d r:%d\n", tagR[r], r);
//             //     // assert(r < right.size());
//             // }
//             // printf("\n");
//             times++;
//             Rs = temp;
//             // printf("Rs:%d\n", Rs);
//             for(auto& x : left[i]){
//                 if(tagR[x] == times - 1){
//                     // update the v_diff and d_diff
//                     for(auto& xl : right[x]) if(xl != i){
//                         ver_intersec[xl]--;
//                         deg_intersec[xl] -= right[x].size();
//                         float s = ver_intersec[xl];
//                         // float s = deg_intersec[xl] * 1.0 / (ver_intersec[xl] + 1);
//                         pq.push(make_pair(s, xl));
//                         v_pq.push(make_pair(ver_intersec[xl], xl));
//                     }
//                 }
//             }
//             while(!v_pq.empty()){
//                 // printf("L");
//                 pair<int,int> t = v_pq.top();
//                 v_pq.pop();
//                 if(tagL[t.second] == 1 && ver_intersec[t.second] == t.first){
//                     delta_max = t.first;
//                     break;
//                 }
//             }
//             threshold = Ls.size() * Rs * 1.0 / (Ls.size() + 1);
//         }
//         for(auto& x : two_hop){
//             tagL[x] = 0;
//         }
//         tagL[i] = 0;
//         vector<int> Rsol;
        
//         for(auto& r : left[i]){
//             if(tagR[r] == times){
//                 Rsol.push_back(r);  
//             }
//         }
//         c += left[i].size() - Rsol.size();
//         for(auto& r : left[i]){
//             tagR[r] = 0;
//         }
//         if(Rsol.size() > 1 && Ls.size() > 1){
//             set<int> left(Ls.begin(), Ls.end());
//             set<int> right(Rsol.begin(), Rsol.end());
//             input.push_back(NODE(mid_count++, left, right));  
//         }
//     }


//     delete[] tagL;
//     delete[] tagR;
//     delete[] ver_intersec;
//     delete[] deg_intersec;

//     printf("%d mid nodes are generated...\n", mid_count);

//     set<pair<int,int>> T;
//     // vector<vector<int>> T(left);
//     for(auto& node : input){
//         for(auto& l : node.left){
//             for(auto& r : node.right){
//                 // T[l][r] = -1;         
//                 T.insert(make_pair(l,r));
//             }
//         }
//     }
    
//     // assert(e -T.size() == c);
//     printf("%lu edges are left... %d \n", e - T.size(), c);
//     // report the compression ratio.
//     // int es = 0;
//     // for(int i = 0, s = left.size();i < s;i++){
//     //     es += left[i].size();  
//     // }
    
//     return input;
// }

// compress the graph directly
void Bi2Tri::duplicate_free_compress(){
    int* tagL = new int[left.size()];
    memset(tagL, 0, sizeof(int) * left.size());

    int* tagR = new int[right.size()];
    memset(tagR, 0, sizeof(int) * right.size());

    // ui* deg_intersec = new ui[left.size()];
    // memset(deg_intersec, 0, sizeof(ui) * left.size());

    ui*  ver_intersec = new ui[left.size()];
    memset(ver_intersec, 0, sizeof(ui) * left.size());

    vector<vector<int> > isEdge(left);
    
    int mid_count = 0;
    // int midsize = 0;
    vector<NODE> input;
    
    int c = 0;
    
    vector<pair<int,int> > L;
    for(int i = 0;i < left.size();i++){
        L.push_back(make_pair(left[i].size(), i));
    }

    // sort(L.begin(), L.end(), compar);

    int* visited = new int[left.size()];
    memset(visited, 0, sizeof(ui) * left.size());

    int* jump = new int[left.size()];
    memset(jump, 0, sizeof(ui) * left.size());

    for(int q = 0, s = L.size();q < s;q++) { 
        int i = L[q].second;
        // printf("%d th ..\n",q);
        // if(jump[i] == 1) continue;
        // extract the two-hop neighbors for each vertex.
        // printf("%d th vertex is processed...\n", i);
        ui degs = 0;
        vector<int> one_hop, two_hop;
        tagL[i] = 1;
        for(int j = 0;j < left[i].size();j++) if(isEdge[i][j] != -1){
            int r = left[i][j];
            tagR[r] = 1;
            for(int k = 0;k < right[r].size();k++) {
                int l = right[r][k];
                if(l != i){
                    tagL[l]++;
                    if(tagL[l] == 2){
                        two_hop.push_back(l);
                        // visited[l] = 1;
                    }
                }
            }
        }
        // printf("%d\n", two_hop.size());

        ui vs = left[i].size();
        priority_queue<pair<int, int>> pq;
        float threshold = 0;

        int delta_max = -1; // neg-infinity
        // calculate the difference for each vertex

        for(auto& x : two_hop){
            int d_count = 0, v_count = 0;
            for(int j = 0;j < left[x].size();j++) if(isEdge[x][j] != -1){
            // for(auto& xr : left[x]) {
                int xr = left[x][j];
                if(tagR[xr] == 1){
                    v_count++;
                    // d_count += right[xr].size();
                }
            }
            // deg_intersec[x] = d_count;
            ver_intersec[x] = v_count;
            delta_max = max(delta_max, v_count);
            pq.push(make_pair(ver_intersec[x], x));   
            // pq.push(make_pair(deg_intersec[x] * 1.0 / (ver_intersec[x] + 1), x));   
            // v_pq.push(make_pair(v_count,x));
        }

        // cerr << "1" << endl;
        // int* start = new int[delta_max + 1];
        // memset(start, -1, sizeof(int) * (delta_max + 1));
        
        // int* next = new int[left.size()];

        // for(auto& x : two_hop){
        //     next[x] = start[ver_intersec[x]];
        //     start[ver_intersec[x]] = x;
        // }

        // printf("%f\n", threshold);
        vector<int> Ls;
        Ls.push_back(i);
        ui Rs = 0;
        int times = 1;
        // printf("3");
        // cerr << "2" << endl;
        // while(delta_max > 0){
        //     // cerr << "4" << endl;
        //     while(delta_max > 0 && start[delta_max] == -1) {
        //         delta_max--;
        //         assert(delta_max >= 0);
        //     }
        //     // cerr << "k" << endl;
        //     if(delta_max <= 0 || delta_max <= Ls.size())
        //         break;
        //     // cerr << "k+1" << endl;
        //     int v = -1;
        //     assert(delta_max >= 0);
        //     // cerr << "5" << endl;
        //     for(v = start[delta_max]; v!= -1;){
        //         int tmp = next[v];
        //         int d = ver_intersec[v];
        //         // assert(d <= left[]);
        //         if(d < delta_max){
        //             // cerr << "6" << endl;
        //             next[v] = start[d];
        //             start[d] = v;
        //         }
        //         else{
        //             // 
        //             // cerr << "7" << endl;
        //             Ls.push_back(v);
        //             assert(v <= left.size());
        //             tagL[v] = 0;
        //             int temp = 0;
        //             for(int j = 0;j < left[v].size();j++) if(isEdge[v][j] != -1){
        //             // for(auto& x : left[p.second]){ 
        //                 int x = left[v][j];
        //                 assert(x < right.size());
        //                 if(tagR[x] == times){
        //                     tagR[x]++;   
        //                     temp++;
        //                 }
        //             }
        //             times++;
        //             Rs = temp;
        //             // vector<int> update; 
        //             for(int j = 0;j < left[i].size();j++) if(isEdge[i][j] != -1){
        //                 int x = left[i][j];
        //                 if(tagR[x] == times - 1){
        //                     for(auto& xl : right[x]) if(xl != i){
        //                         ver_intersec[xl]--;
        //                     }
        //                 }
        //             }
        //             // threshold = Ls.size() * Rs * 1.0 / (Ls.size() + 1);
        //             start[delta_max] = tmp;
        //             break;
        //         }
        //         v = tmp;
        //     }
        //     if(v == -1) start[delta_max] = -1;
        // }
        // cerr << "q" << endl;
        // delete[] start;
        // delete[] next;
        // cerr << "q+1" << endl;
        while(!pq.empty()){
            // printf("%d\n", pq.size());
             // update the thresho ld.
            pair<int, int>  p = pq.top();
            pq.pop();
            if(p.first <= (int)floor(threshold)){
                break;
            }
            // pair<float, int>  p = pq.top();
            // pq.pop();
            // assert(p.second < left.size());
            // if(fabsf(p.first - ver_intersec[p.second]) > EPSINON){
            //     continue;
            // }
            if(p.first != ver_intersec[p.second]){
                continue;
            }
            // if(fabsf(p.first - deg_intersec[p.second] * 1.0 / (ver_intersec[p.second] + 1)) > EPSINON){
            //     continue;
            // }
            // printf("merge %d\n", p.second);
            Ls.push_back(p.second);
            tagL[p.second] = 0;
            int temp = 0;
            for(int j = 0;j < left[p.second].size();j++) if(isEdge[p.second][j] != -1){
            // for(auto& x : left[p.second]){ 
                int x = left[p.second][j];
                if(tagR[x] == times){
                    tagR[x]++;   
                    temp++;
                }
            }
            times++;
            Rs = temp;
            // printf("Rs:%d\n", Rs);
            
            vector<int> update; 
            for(int j = 0;j < left[i].size();j++) if(isEdge[i][j] != -1){
            // for(auto& x : left[i]){
                int x = left[i][j];
                if(tagR[x] == times - 1){
                    // update the v_diff and d_diff
                    for(auto& xl : right[x]) if(xl != i){
                        ver_intersec[xl]--;
                        // deg_intersec[xl] -= right[x].size();
                        if(visited[xl] == 0){
                            update.push_back(xl);
                            visited[xl] = 1;
                        }
                        // float s = ver_intersec[xl];
                        // pq.push(make_pair(s, xl));
                        // v_pq.push(make_pair(ver_intersec[xl], xl));
                    }
                }
            }

            for(auto& xl : update){
                int s = ver_intersec[xl];
                pq.push(make_pair(s, xl));
                // v_pq.push(make_pair(ver_intersec[xl], xl));
                visited[xl] = 0;
            }

            
            // while(!v_pq.empty()){
            //     // printf("L");
            //     pair<int,int> t = v_pq.top();
            //     v_pq.pop();
            //     if(tagL[t.second] != 0  && ver_intersec[t.second] == t.first){
            //         delta_max = t.first;
            //         break;
            //     }
            // }
            threshold = Ls.size() * Rs * 1.0 / (Ls.size() + 1);
            // threshold = Ls.size() * 1.0;
        }

        // cerr << "3" << endl;
        
        for(auto& x : two_hop){
            tagL[x] = 0;
            // visited[x] = 0;
        }
        tagL[i] = 0;
        // visited[i] = 0;
        vector<int> Rsol, Ridx;
        
        for(int j = 0;j < left[i].size();j++){
            int r = left[i][j];
            if(tagR[r] == times){
                Rsol.push_back(r);  
            }
        }

        // c += left[i].size() - Rsol.size();
        
        if(Rsol.size() > 1 && Ls.size() > 1){
            // printf("%d %dedges are compressed\n", Rsol.size() , Ls.size());
            set<int> left_(Ls.begin(), Ls.end());
            set<int> right_(Rsol.begin(), Rsol.end());
            // printf("%d %d\n",right_.size() , left_.size());
            input.push_back(NODE(mid_count++, left_, right_));
            for(auto& ls : Ls){
                jump[ls] = 1;
                for(int j = 0;j < left[ls].size();j++){
                    int r = left[ls][j];
                    if(tagR[r] == times && isEdge[ls][j] != -1){
                        // printf("delete: %d %d\n", ls, r);
                        isEdge[ls][j] = -1;
                    }
                }
            }  
        }

        for(auto& r : left[i]){
            tagR[r] = 0;
        }
    }

    delete[] tagL;
    delete[] tagR;
    delete[] ver_intersec;
    delete[] jump;
    // delete[] deg_intersec;
    delete[] visited;
    

    printf("%d mid nodes are generated...\n", mid_count);

    int ec = 0;
    for(auto& node : input){
        ec += node.left.size() + node.right.size();
        // printf("L:");
        // for(auto& l : node.left){
        //     printf("%d ", l);
        // }
        // printf("R:");
        // for(auto& r : node.right){
        //     printf("%d ", r);
        // }
        // printf("\n");
    }
    // printf("ec:%d\n", ec);

    for(int j = 0;j < isEdge.size();j++){
        auto& l = isEdge[j];
        for(auto& r : l){
            if(r != -1){
                // printf("%d %d\n", j, r);
                ec++;
            }
        }
    }
    
    printf("compression ratio: %f%\n", ec * 100.0 / this->e);
    // assert(e -T.size() == c);
    // printf("%lu edges are left... %d \n", e - T.size(), c);
    // report the compression ratio.
    // int es = 0;
    // for(int i = 0, s = left.size();i < s;i++){
    //     es += left[i].size();  
    // }
    
}


void Bi2Tri::summary(){
    // give the summary of current graph. 
    int L = 0, R = 0;
    for(int i = 0;i < ln;i++){
        L += L_valid[i];
    }

    for(int i = 0;i < rn;i++){
        R += R_valid[i];
    }

    // for(auto l : L_valid){
    //     L += l;
    // }
    // for(auto r : R_valid){
    //     R += r;
    // }
    int e = 0;
    for(int i = 0,s = left.size();i < s;i++){
        if(L_valid[i] == 1){
            for(auto n : left[i]){
                if(R_valid[n] == 1){
                    e++;
                }
            }
        }
    }
    printf("|L| : %d, |R| : %d, |E| : %d \n", L, R, e);
}

