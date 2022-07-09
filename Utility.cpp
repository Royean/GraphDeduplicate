#include "Utility.h"

int Rand(int i){
    return rand()%i;
}

Bucket_queue::Bucket_queue(int max_priority, int num_of_elements){
    bucket = new int[max_priority + 1];
    next = new int[num_of_elements];
    memset(bucket, -1, sizeof(int) * (max_priority + 1));
    memset(next, -1, sizeof(int) * (num_of_elements));
}


// Bucket_queue::~Bucket_queue(){
//     if(this->next != nullptr){
//         delete[] next;
//     }
//     if(this->bucket != nullptr){
//         delete[] bucket;
//     }
// }

void Bucket_queue::build_queue(vector<pair<int,int>>& data){
    for(auto& p : data){
        next[p.second] = bucket[p.first];
        bucket[p.second] = p.first;
    }
}

bool Bucket_queue::is_empty(){
    return (this->current_priority == 0) && (this->bucket[current_priority] == -1);
}

void Bucket_queue::adjust(int new_priority, int element_id){
   next[element_id] = bucket[new_priority];
   bucket[new_priority] = element_id;
}

int Bucket_queue::pop(){
    // return bucket[current_priority];
    // // while(current_priority > 0){
    //     // while(current_priority >= 0 && bucket[current_priority] == -1) current_priority--;
    //     if(current_priority < 0){
    //         // break;
    //         return -1;
    //     }
    //     while(bucket[current_priority] == -1){
    //         current_priority--;
    //     }
    //     int ans = bucket[current_priority];
    //     int ans =  bucket[current_priority];

    //     while(current_priority >= 0 && bucket[current_priority] == -1){
    //         current_priority--;

    //     } 
    //     int v = -1;
    //     for(v = bucket[current_priority]; v != -1;){
    //         int tmp = next[v];

    //     }
    // }

    // while (max_d >= 0 && bin_head[max_d] == -1) max_d--;
	// 		if (max_d < 0) break;
	// 		int v = -1;
	// 		for (v = bin_head[max_d]; v != -1;) {
	// 			int tmp = bin_next[v];
	// 			int deg = head[v].neiSum - head[v].weight;
	// 			if (head[v].del == 0 && deg > 0) {
	// 				if (deg < max_d) {
	// 					bin_next[v] = bin_head[deg];
	// 					bin_head[deg] = v;
	// 				}
	// 				else {
	// 					deleteVertex(v, pone, degree_one, degree_two);
	// 					backTrack.push_back(make_pair(v, n));
	// 					bin_head[max_d] = tmp;
	// 					break;
	// 				}
	// 			}
	// 			v = tmp;
	// 		}
	// 		if (v == -1) bin_head[max_d] = -1;

}

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

vector<string> split_line(string line, string delimiter = "\r\t\n "){
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