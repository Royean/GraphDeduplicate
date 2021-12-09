//
//  FPGrowth.cpp
//  FrequentPatternMining


#include "FPGrowth.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#pragma -mark public method

FPGrowth::FPGrowth(float _minSupport)
{
    //init root node
    root = new FPTreeNode("",NULL,0);
    //set min support
    minSupport = _minSupport;
}

void FPGrowth::initFromFile(string fileName)
{
    ifstream ifs("./bigraph/"+ fileName, std::ifstream::in);
    int e = 0;
    if (ifs) {
        string line;
        getline(ifs, line);
        left.resize(atoi(line.c_str()));
        getline(ifs, line);
        right.resize(atoi(line.c_str()));
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
                int temp = atoi(content[k++].c_str());
                left[l].push_back(temp);
                right[temp].push_back(l);
                e++;
            }
        }
    }
    ifs.close();
    this->e = e;

    cout << "Init from file: " << fileName  << " ,with min support: " << minSupport <<endl;
    flush(cout);
    start = clock();
    buildHeaderTable(fileName);
    buildFPTree(fileName);
    prefix = NULL; //insure that first mining will creat a new prefix
}

bool ResultComp(const FPFreqResult& a, const FPFreqResult& b){
    return (a.freq - 1) * (a.items.size() - 1) > (b.freq - 1) * (b.items.size() - 1);
}

void FPGrowth::outputToFile(string fileName)
{
    ofstream ofile;
    ofile.open(fileName.c_str());
    if (!ofile) {
        cout << "Can not open file :" << fileName << endl;
        return;
    }
    
    // vector<int> 

    // vector<FPHeaderTableNode *> FPHeaderTable;
    vector<FPHeaderTableNode *> FPHeaderTable_Copy(FPHeaderTable);
    printf("header size: %d\n", FPHeaderTable_Copy.size());

    // int cluster = 3;

    int inc = FPHeaderTable_Copy.size();
    
    set<pair<int,int>> ET;
    int ec = 0;
    for(int k = 0;k < 1;k++){
        vector<FPFreqResult> result;
        vector<FPHeaderTableNode *> temp;
        for(int s = k * inc, j = s; j < FPHeaderTable_Copy.size() && j < s + inc;j++){
            temp.push_back(FPHeaderTable_Copy[j]);
        }
        // printf("temp size: %d\n", temp.size());
        FPHeaderTable = temp;
        output(result);
        int vc = 0;
        int m = 0;
        // printf("result size: %d\n", result.size());
        // sort
        sort(result.begin(), result.end(), ResultComp);
        
        for(auto& it : result){
           
            if(it.items.size() > 1 && it.freq > 2) {
                //  printf(" %d / %d\n",it.items.size(), it.freq);
                // cout << "success" << endl;
                // check if exist duplicate path.
                vc++;
                vector<FPTransItem>::iterator itemIt;
                vector<int> R;
                for (itemIt = it.items.begin(); itemIt != it.items.end(); itemIt ++) {
                    // cout << itemIt->name << endl;
                    R.push_back(atoi(itemIt->name.c_str()));
                    // cout << R.back() << endl;
                }

                // Res represents left vertices.
                vector<int> Res;
                int* tag = new int[left.size()];
                memset(tag, 0, sizeof(int) * left.size());

                vector<int> Temp;
                for(int k = 0;k < R.size();k++){
                    int r = R[k];
                    for(auto& l : right[r]){
                        tag[l]++;
                        Temp.push_back(l);
                        // printf("%d %d\n", l, tag[l]);
                        if(tag[l] == it.items.size()){
                            Res.push_back(l);
                        }
                    }
                }

                for(auto& x : Temp){
                    tag[x] = 0;
                }

                bool ok = true;
                for(auto& l : Res){
                    for(auto& r : R){
                        pair<int, int> p = make_pair(l,r);
                        if(ET.find(p) != ET.end()){
                            ok = false;
                        }
                    }
                }
                if(ok){
                    // printf("L: ");
                    // for(auto& l : Res){
                    //     printf("%d ", l);
                    //     // ET.insert(make_pair(l,r));
                    // }
                    // printf("R: ");
                    // for(auto& r : R){
                    //     printf("%d ", r);
                    //     // ET.insert(make_pair(l,r));
                    // }
                    // printf("\n");
                    for(auto& l : Res){
                        for(auto& r : R){
                            pair<int, int> p = make_pair(l,r);
                            // if(ET.find(p) != ET.end()){
                            //     ok = false;
                            // }
                            ET.insert(make_pair(l,r));
                        }
                    }
                    ec += Res.size() + R.size();
                }
            }
        }
        
    }

    printf("Compression ratio : %f\n", (e - ET.size() + ec) * 1.0 / e);
    time_t t_end = time(NULL);
    cout<< "Total time used for mining: "<< (clock()-start)/double(CLOCKS_PER_SEC) << " (sec)" << endl;
    flush(cout);
    vector<FPFreqResult>::iterator it;

    

    // for (it = result.begin(); it != result.end(); it ++) {
    //     vector<FPTransItem>::iterator itemIt;
    //     for (itemIt = it->items.begin(); itemIt != it->items.end(); itemIt ++) {
    //         ofile << itemIt->name << " " ;
    //     }
    //     ofile << " : " << it->freq << endl;
    // }
    // t_end = time(NULL);
    cout<< "Total time used for mining and output: "<< (clock()-start)/double(CLOCKS_PER_SEC) << " (sec)" << endl;
    flush(cout);
    ofile.close();
}

void FPGrowth::output(vector<FPFreqResult> &result)
{
    miningFPTree(result);
}

//any allocated space in trans need caller to delete immediately
void FPGrowth::initFromTrans(vector<FPTrans> &trans)
{
    buildHeaderTable(trans);
    buildFPTree(trans);
}

#pragma -mark private method

extern vector<string> split_line(string line, string delimiter);

void FPGrowth::buildHeaderTable(string fileName)
{
    ifstream ifile;
    ifile.open(("./bigraph/"+ fileName).c_str());
    if (!ifile) {
        cout << "Can not open file :" << fileName << endl;
        return;
    }
    
    string tran; //get a transaction
    int index = 0; //index in HeaderTable
    int transCount = 0;
    getline(ifile,tran);
    getline(ifile,tran);
    while (getline(ifile,tran)) {
        transCount ++;
        // cout << tran << endl;
        stringstream strStream(tran);
        string item;
        strStream >> item;
        strStream >> item;
        while (strStream >> item) {
            if (nameIndex.find(item) != nameIndex.end()) {
                FPHeaderTable[nameIndex[item]]->freq ++;
            }
            else{ //item doesn't exit in Header Table
                nameIndex[item] = index;
                FPHeaderTableNode *headerTableNode = new FPHeaderTableNode(item);
                FPHeaderTable.push_back(headerTableNode);
                index ++;
            }
        }
    }
    ifile.close();
    
    
    //delete items that not match min Support
    vector<FPHeaderTableNode *>::iterator it;
    // minCountSupport = (float)transCount * minSupport;
    minCountSupport = 2.1f;
    for (it = FPHeaderTable.begin(); it != FPHeaderTable.end(); it++)
    {
        FPHeaderTableNode *node = *it;
        if ((float)(node->freq) < minCountSupport) {
            it = FPHeaderTable.erase(it);
            it --;
        }
    }
    
    sortHeaderTable();
}

void FPGrowth::buildHeaderTable(vector<FPTrans> &trans)
{
    int index = 0; //index in HeaderTable
    int transCount = 0;
    vector<FPTrans>::iterator transIt;
    for (transIt = trans.begin(); transIt != trans.end(); transIt++) {
        transCount ++;
        vector<FPTransItem>::iterator itemIt;
        for (itemIt = transIt->items.begin(); itemIt != transIt->items.end(); itemIt++) {
            if (nameIndex.find((itemIt)->name) != nameIndex.end()) {
                FPHeaderTable[nameIndex[(itemIt)->name]]->freq += (itemIt)->count;
            }
            else{ //item doesn't exit in Header Table
                nameIndex[(itemIt)->name] = index;
                FPHeaderTableNode *headerTableNode = new FPHeaderTableNode((itemIt)->name,(itemIt)->count);
                FPHeaderTable.push_back(headerTableNode);
                index ++;
            }
        }
    }
    
    //delete items that not match min Support
    vector<FPHeaderTableNode *>::iterator it;
    
    for (it = FPHeaderTable.begin(); it != FPHeaderTable.end(); it++)
    {
        FPHeaderTableNode *node = *it;
        if ((float)(node->freq) < minCountSupport) {
            it = FPHeaderTable.erase(it);
            it --;
        }
    }
    
    sortHeaderTable();
}


bool CompareHeaderTableNode(FPHeaderTableNode *a,FPHeaderTableNode *b)
{
    return a->freq > b->freq;
}

void FPGrowth::sortHeaderTable()
{
    sort(FPHeaderTable.begin(), FPHeaderTable.end(), CompareHeaderTableNode);

    //update name index
    int index = 0;
    nameIndex.clear();
    vector<FPHeaderTableNode *>::iterator it;
    for (it = FPHeaderTable.begin(); it != FPHeaderTable.end(); it++)
    {
        FPHeaderTableNode *node = *it;
        nameIndex[node->name] = index;
        index ++;
    }
    itemCount = index;
}


void FPGrowth::insertTran(vector<string> &items)
{
    FPTreeNode *currentTreeNode = root;
    unordered_map<string, FPTreeNode*> *currentChildren;
    currentChildren = &(currentTreeNode->children);
    vector<string>::iterator it;
    for (it = items.begin(); it != items.end(); it ++) {
        
        //if prefix match
        if (currentChildren->find(*it) != currentChildren->end()) {
            (*currentChildren)[*it]->count ++;
            currentTreeNode = (*currentChildren)[*it];
            currentChildren = &(currentTreeNode->children);
        }
        else{
            FPTreeNode *newTreeNode = new FPTreeNode(*it,currentTreeNode);
            (*currentChildren)[*it] = newTreeNode;
            newTreeNode->parent = currentTreeNode; //set parent
            currentTreeNode = newTreeNode;
            currentChildren = &(newTreeNode->children);
            
            //update FPHeaderTable
            if (FPHeaderTable[nameIndex[*it]]->head == NULL) {
                FPHeaderTable[nameIndex[*it]]->head = currentTreeNode;
                FPHeaderTable[nameIndex[*it]]->end = currentTreeNode;
            }
            else {
                FPHeaderTable[nameIndex[*it]]->end->link = currentTreeNode;
                FPHeaderTable[nameIndex[*it]]->end = currentTreeNode;
            }
        }
        

    }
}

void FPGrowth::insertTran(vector<FPTransItem*> &items)
{
    FPTreeNode *currentTreeNode = root;
    unordered_map<string, FPTreeNode*> *currentChildren;
    currentChildren = &(currentTreeNode->children); //a map from item name to children tree node
    
    vector<FPTransItem*>::iterator it;
    for (it = items.begin(); it != items.end(); it ++) {
        
        //if prefix match, search next
        if (currentChildren->find((*it)->name) != currentChildren->end()) {
            
            (*currentChildren)[(*it)->name]->count += (*it)->count;
            
            currentTreeNode = (*currentChildren)[(*it)->name];
            currentChildren = &(currentTreeNode->children);
        }
        else{ //prefix not match
            FPTreeNode *newTreeNode = new FPTreeNode((*it)->name,currentTreeNode);
            
            newTreeNode->count = (*it)->count;
            
            (*currentChildren)[(*it)->name] = newTreeNode;
            newTreeNode->parent = currentTreeNode; //set parent
            currentTreeNode = newTreeNode;
            currentChildren = &(newTreeNode->children);
            
            //update FPHeaderTable
            if (FPHeaderTable[nameIndex[(*it)->name]]->head == NULL) {
                FPHeaderTable[nameIndex[(*it)->name]]->head    = currentTreeNode;
                FPHeaderTable[nameIndex[(*it)->name]]->end     = currentTreeNode;
            }
            else {
                FPHeaderTable[nameIndex[(*it)->name]]->end->link   = currentTreeNode;
                FPHeaderTable[nameIndex[(*it)->name]]->end         = currentTreeNode;
            }
        }

    }
}

//sort a transaction use nameIndex
unordered_map<string, int> *gNameIndex = NULL;
bool CompareItem(string a,string b)
{
    return (*gNameIndex)[a] < (*gNameIndex)[b];
}
bool CompareTransItem(FPTransItem *a,FPTransItem *b)
{
    return (*gNameIndex)[a->name] < (*gNameIndex)[b->name];
}

void FPGrowth::buildFPTree(string fileName)
{
    ifstream ifile;
    ifile.open(("./bigraph/"+ fileName).c_str());
    if (!ifile) {
        cout << "Can not open file :" << fileName << endl;
        return;
    }
    string tran; //get a transaction
    int index = 0; //index in HeaderTable
    int transCount = 0;
    getline(ifile,tran);
    getline(ifile,tran);
    while (getline(ifile,tran)) {
        transCount ++;
        stringstream strStream(tran);
        string item;
        vector<string> items; //transaction items
        strStream >> item;
        strStream >> item;
        while (strStream >> item) {
            if (nameIndex.find(item) != nameIndex.end()) {
                items.push_back(item); //only add item that match minSupport
            }
        }
        
        //sort items
        gNameIndex = &nameIndex;
        sort(items.begin(), items.end(), CompareItem);
        gNameIndex = NULL;
        
        insertTran(items);//insert a modified transaction to FPTree
    }
    ifile.close();
}

void FPGrowth::buildFPTree(vector<FPTrans> &trans)
{
    vector<FPTrans>::iterator transIt;
    vector<FPTransItem>::iterator itemIt;
    for (transIt = trans.begin(); transIt != trans.end(); transIt++) {
        vector<FPTransItem*> items;
        for (itemIt = transIt->items.begin(); itemIt != transIt->items.end(); itemIt++) {
            if (nameIndex.find((itemIt)->name) != nameIndex.end()) {
                items.push_back(&(*itemIt)); //only add item that match minSupport
            }
        }
        //sort items
        gNameIndex = &nameIndex;
        sort(items.begin(), items.end(), CompareTransItem);
        gNameIndex = NULL;
        insertTran(items);//insert a modified transaction to FPTree
    }
}

#pragma -mark  Mining
void FPGrowth::printSingleResult(vector<FPFreqResult> &result)
{
    FPFreqResult resultLine;
    vector<FPTransItem *>::iterator it;
    it = prefix->end() - 1;
    int supportCount = (*it)->count;
    for (; it >= prefix->begin(); it --) {
        FPTransItem item((*it)->name,0);
        resultLine.items.push_back(item);
    }
    resultLine.freq = supportCount;
    if(resultLine.items.size() > 1 && supportCount > 1){
        result.push_back(resultLine);
    }
}

void FPGrowth::miningFPTree(vector<FPFreqResult> &result)
{
    // 
    if (prefix == NULL){ //First mining
        prefix = new vector<FPTransItem *>;
    }
    
    // && result.size() > (int)(left.size() * 0.7 / 1)
    if (FPHeaderTable.size() == 0) { //mining arrives root
        if (prefix->size() > 0) {
            printSingleResult(result);
        }
        if(result.size() >= (int)(0.7 * left.size()))
        return;
    }
    
    if (prefix->size() > 0) { //print current prefix
        printSingleResult(result);
    }
    
    //mining each item in HeaderTable
    vector<FPHeaderTableNode *>::iterator it;
    
    for (it = FPHeaderTable.end() - 1; it >= FPHeaderTable.begin(); it--)
    {
        FPTransItem *newPrefix = new FPTransItem((*it)->name,0);
        vector<FPTrans> trans;
        FPTreeNode *cur;
        for (cur = (*it)->head; cur != (*it)->end; cur = cur->link) {
            vector<FPTransItem> items;
            //search from node to root,build a modified trans
            FPTreeNode *curToRoot;
            int basicCount = cur->count;
            for (curToRoot = cur->parent;curToRoot->parent != NULL;curToRoot = curToRoot->parent){
                FPTransItem item(curToRoot->name,basicCount);
                items.push_back(item);
            }
            FPTrans tranLine = FPTrans(items);
            trans.push_back(tranLine);
        }
        
        //when cur == (*it)->end
        vector<FPTransItem> items;
        //search from node to root,build a modified trans
        FPTreeNode *curToRoot;
        int basicCount = cur->count;
        for (curToRoot = cur->parent;curToRoot->parent != NULL;curToRoot = curToRoot->parent){
            FPTransItem item(curToRoot->name,basicCount);
            items.push_back(item);
        }
        FPTrans tranLine = FPTrans(items);
        trans.push_back(tranLine);
        
        newPrefix->count = (*it)->freq;
        
        //add prefix
        prefix->push_back(newPrefix);
        
        //recursive mining
        FPGrowth *fpTemp = new FPGrowth(minSupport);
        fpTemp->minCountSupport = 2;
        fpTemp->initFromTrans(trans);
        fpTemp->prefix = prefix;

        fpTemp->miningFPTree(result);
        
        delete fpTemp;
        prefix->pop_back();
    }

}

#pragma -mark destructor
FPGrowth::~FPGrowth()
{
    if (root != NULL) {
        delete root;
    }
    if (prefix != NULL && prefix->size() == 0) {
        vector<FPTransItem *>::iterator it;
        for (it = prefix->begin(); it != prefix->end(); it++)
        {
            if (*it!=NULL) {
                delete *it;
            }
        }
        delete prefix;
    }
    vector<FPHeaderTableNode *>::iterator it;
    for (it = FPHeaderTable.begin(); it != FPHeaderTable.end(); it++)
    {
        delete *it;
    }
}


#pragma -mark debug method

void FPGrowth::outputHeaderTable()
{
    cout << "HeaderTable: "<< endl;
    vector<FPHeaderTableNode *>::iterator it;
    for (it = FPHeaderTable.begin(); it != FPHeaderTable.end(); it++)
    {
        FPHeaderTableNode *node = *it;
        cout << "name:" << node->name
            << "\tfreq: "<< node->freq
            << "\thead: "<<node->head
            << "\tend: "<<node->end
            <<endl;
    }
    cout << endl;
}

void FPGrowth::outputTran(vector<string> items)
{
    vector<string>::iterator it;
    cout << "Tran items: " ;
    for (it = items.begin(); it != items.end(); it ++) {
        cout << *it << " ";
    }
    cout << endl;
}

void FPGrowth::outputTree()
{
    FPTreeNode* cur = root;
    queue<FPTreeNode *>q;
    q.push(cur);
    cout << "Current Tree is:" << endl;
    while (!q.empty()) {
        cur = q.front();
        q.pop();
        cout << "Name: " << cur->name << "\tCount: " << cur->count << "\tChildren num: "<< cur->children.size() << endl;
        unordered_map<string, FPTreeNode*>::iterator it;
        for (it = cur->children.begin(); it != cur->children.end(); it ++) {
            q.push((*it).second);
        }
    }
    cout << endl;
}

void FPGrowth::outputPrefix()
{
    cout << "Current Prefix is: "<<endl;
    vector<FPTransItem *>::iterator it;
    it = prefix->end() - 1;
    int supportCount = (*it)->count;
    for (; it >= prefix->begin(); it --) {
        cout << (*it)->name << " ";
    }
    cout << supportCount << endl;
}

void FPGrowth::outputFreq(){
    vector<FPTransItem *>::iterator it;
    it = prefix->end() - 1;
    int supportCount = (*it)->count;
    cout << "================================" << endl;
    for (; it >= prefix->begin(); it --) {
        cout << (*it)->name << " ";
    }
    cout << supportCount << endl;
    cout << "================================" << endl << endl;
}
