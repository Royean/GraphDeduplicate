#include "Utility.h"

int Rand(int i){
    return rand()%i;
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