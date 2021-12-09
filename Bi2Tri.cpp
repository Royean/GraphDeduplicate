#pragma once
#include "Bi2Tri.h"

bool compar(const pair<int,int>& a, const pair<int,int>& b){
    if(a.first < b.first) return true;
    return false;
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

Bi2Tri::Bi2Tri(int n, int m, float prob){
    // int n = n, m;
    left.resize(n);
    right.resize(m);

    vector<int> L;
    
    for(int i = 0;i < m;i++){
        L.push_back(i);
    }

    // int e = 0;
    srand(unsigned(time(0)));
    for(int i = 0;i < n;i++){
        random_shuffle(L.begin(), L.end(), Rand);
        for(int j = 0;j < (int)(prob * m);j++){
            left[i].push_back(L[j]);
            right[L[j]].push_back(i);
            e++;
        }
    }
    printf("bigraph reading done......\n");
    printf("left nodes: %ld, right nodes: %ld, edges: %d\n", left.size(), right.size(), e);
}


Bi2Tri::Bi2Tri(string filename){
    ifstream ifs("./bigraph/"+ filename, std::ifstream::in);
    int e = 0;
    if (ifs) {
        string line;
        getline(ifs, line);
        left.resize(atoi(line.c_str()));
        getline(ifs, line);
        right.resize(atoi(line.c_str()));
        printf("1\n");
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
    printf("bigraph reading done......\n");
    printf("left nodes: %ld, right nodes: %ld, edges: %d\n", left.size(), right.size(), e);
}

void Bi2Tri::reduce(){
    vector<int> L_degree_one, R_degree_one;
    for(int i = 0, s = left.size();i < s;i++){
        if(left[i].size() == 1){
            L_degree_one.push_back(i);
        }
    }
    
    for(int i = 0, s = left.size();i < s;i++){
        if(right[i].size() == 1){
            R_degree_one.push_back(i);
        }
    }

    while(!L_degree_one.empty() || !R_degree_one.empty()){
        while (!L_degree_one.empty()){
            

        }
        while (!R_degree_one.empty()){
            
        }
    }
}

vector<NODE> Bi2Tri::convert(int redundancy){
    int* tagL = new int[left.size()];
    memset(tagL, 0, sizeof(int) * left.size());

    int* tagR = new int[right.size()];
    memset(tagR, 0, sizeof(int) * right.size());

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

    sort(L.begin(), L.end(), compar);

    int* visited = new int[left.size()];
    memset(visited, 0, sizeof(ui) * left.size());

    int* jump = new int[left.size()];
    memset(jump, 0, sizeof(ui) * left.size());

    for(int q = 0, s = L.size();q < s;q++){ 
        int i = L[q].second;
        // if(jump[i] == 1){
        //     continue;
        // } 
        // extract the two-hop neighbors for each vertex.
        // printf("%d th vertex is processed...\n", i);
        ui degs = 0;
        vector<int> one_hop, two_hop;
        tagL[i] = 1;
        for(int j = 0;j < left[i].size();j++) {
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
            // printf("Rs:%d\n", Rs);
            
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
            if(redundancy > 1){
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

                redundancy--;
            }

            threshold = Ls.size() * Rs * 1.0 / (Ls.size() + 1);
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
void Bi2Tri::compress(){
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
    // printf("left:%d right:%d\n", left.size(), right.size());
    // for(auto& r: left){
    //     for(auto& e : r){
    //         printf("%d ", e);
    //     }
    //     printf("\n");
    // }
}

