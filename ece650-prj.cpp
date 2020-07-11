// Compile with c++ ece650-a2cpp -std=c++11 -o ece650-a2
// note: add upper time limit
//
#include <iostream>
#include <memory>
#include <string>
#include <regex>
#include <stdexcept>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <pthread.h>
#include <semaphore.h>
#include <ctime>
#include <iomanip>
#include <chrono>
#include<mutex>
#include <sstream>
//defined std::unique_ptr
#include <memory>
// defines Var and Lit
#include "minisat/core/SolverTypes.h"
// defines Solver
#include "minisat/core/Solver.h"

#include "a2.h"

using namespace Minisat;
using namespace std;
using namespace chrono;
const int INF = 10000;
typedef pair<int,int> int_pairs;

typedef struct thread_data{
    Edge_list* thread_edge_list;
    vector<int> cnf_result;
    vector<int> approx_vc_1_result;
    vector<int> approx_vc_2_result;
    int cnf_time = 0;
    int vc_1_time = 0;
    int vc_2_time = 0;

} thread_data;

const static int MAX_TIME = 2000000000; ///set terminalte time: 30 sec;
time_t start_time, run_time;
sem_t sem0, sem2, sem3, sem1, sem4;//信号量
time_point<high_resolution_clock> start1;
time_point<high_resolution_clock> end1;
mutex p_mutex;

/*to update:
 * 1, reduce variable number to the correct index
 * for example, enter V 35 and E {<35,34>,<34,22>,<22,1>}. There is actually three vertex instead of 35
 * create a array X to store the correct index and run CNF using the index of X.
*/

vector<int> getCNF(const Edge_list& temp_list){
    unique_ptr<Solver> solver(new Solver());
    int n = temp_list.get_v();// # of vertex
    int k = 1;
    vector<int> result; //store result
    map<int_pairs,Lit> my_map; // use this map to store used x, there are a total of n*k possible x clause
    // pair format <n,k>
    map<int_pairs,Lit>::iterator itr_p;
    map<int_pairs,Lit>::iterator itr_q;

//    // 优化
//    Edge_list new_temp_list = temp_list;
//    vector<int> temp_vec;
//    vector <int> unique_vertex_list;
//    Edge temp_edge2;
//
//    while(!new_temp_list.empty()){
//        temp_edge2 = new_temp_list.pop();
//        temp_vec.push_back(temp_edge2.start);
//        temp_vec.push_back(temp_edge2.end);
//    }
//
//    set<int> temp_set(temp_vec.begin(), temp_vec.end());
//    unique_vertex_list.assign(temp_set.begin(),temp_set.end());
//
//    new_temp_list = temp_list;
//    if(unique_vertex_list.size() < temp_list.get_v()){//replace edge with reduced number
//        n = unique_vertex_list.size();
//        for (int i = 0; i < temp_list.size(); i++){
//            temp_edge2 = new_temp_list.pop();
//            for (int j = 0; j < unique_vertex_list.size(); j++){
//                if (temp_edge2.start == unique_vertex_list[j]){
//                    temp_edge2.start = j;
//                }
//                else if(temp_edge2.end == unique_vertex_list[j]){
//                    temp_edge2.end = j;
//                }
//            }
//            new_temp_list.push(temp_edge2);
//        }
//    }
//    // 结束

    start_time = clock();//set start time;

    while (k<n){ //k = [1...n-1]. the vertex cover size
        run_time = clock();
        if(run_time - start_time > MAX_TIME){
            result.clear();
            break;
        }
        // clause 1:
        for (int counter_k = 1; counter_k <= k; counter_k ++){ //∀i∈[1,k]
            //∀i ∈ [1, k], a clause (x 1,i ∨ x 2,i ∨ · · · ∨ x n,i )
            vec<Lit> my_clause1;//append all literal below to this vector.
            for (int counter_n = 0; counter_n < n; counter_n ++){
                int_pairs temp_key = make_pair(counter_n, counter_k);
                Lit temp_lit_1{};
                temp_lit_1 = mkLit(solver->newVar());
                my_map.insert(pair<int_pairs,Lit>(temp_key,temp_lit_1));
                my_clause1.push(temp_lit_1);
            }
            solver->addClause_(my_clause1);
        }
        //clause 2:
        if (k > 1){
            for (int counter_m = 0; counter_m < n; counter_m ++){
                for (int counter_p = 1; counter_p < k; counter_p++){// p [1...k-1]
                    int_pairs temp_p = make_pair(counter_m, counter_p);
                    itr_p = my_map.find(temp_p);

                    for (int counter_q = counter_p+1; counter_q <= k; counter_q++) {// q [p...k]
                        int_pairs temp_q = make_pair(counter_m, counter_q);
                        itr_q = my_map.find(temp_q);
                        solver->addClause(~itr_p->second,~itr_q->second);
                    }
                }
            }
        }
        //clause 3:
        if (n > 1){
            for (int counter_m = 1; counter_m <= k; counter_m ++){
                for (int counter_p = 0; counter_p < n-1; counter_p++){ // p [1...n-1]
                    int_pairs temp_p = make_pair(counter_p, counter_m);
                    itr_p = my_map.find(temp_p);
                    for (int counter_q = counter_p+1; counter_q < n; counter_q ++){ // q [p+1 ... n]
                        int_pairs temp_q = make_pair(counter_q, counter_m);
                        itr_q = my_map.find(temp_q);
                        solver->addClause(~itr_p->second,~itr_q->second);
                    }
                }
            }
        }

        //clause 4:
        Edge_list my_list = temp_list;
//        ∀<i, j> ∈ E, a clause (x i,1 ∨ x i,2 ∨ · · · ∨ x i,k ∨ x j,1 ∨ x j,2 ∨ · · · ∨ x j,k )
        while(!my_list.empty()){ // loop through edge list to ensure every edge is incident to at least one vertex in the vertex cover
            Edge temp_edge = my_list.pop();//            cout<<temp_edge.start<<", "<<temp_edge.end<<endl;
            vec<Lit> my_clause4;//append all literal below to this vector.
            for (int counter_k = 1; counter_k <= k; counter_k++){
                int_pairs temp_i = make_pair(temp_edge.start, counter_k);
                itr_p = my_map.find(temp_i);
                my_clause4.push(itr_p->second);
            }
            for (int counter_k = 1; counter_k <= k; counter_k++){
                int_pairs temp_j = make_pair(temp_edge.end, counter_k);
                itr_q = my_map.find(temp_j);
                my_clause4.push(itr_q->second);
            }
            solver->addClause_(my_clause4);
        }
//        my_list.~Edge_list();


        // clause 5
        for(int counter_i = 0; counter_i < n; counter_i++){
            for (int counter_j = 1; counter_j < k; counter_j ++){
                for (int i = 0; i <counter_i; i ++){
                    int_pairs temp_p = make_pair(counter_i, counter_j);
                    itr_p = my_map.find(temp_p);
                    for (int j = counter_j; j <= k; j++){
                        int_pairs temp_q = make_pair(i, j);
                        itr_q = my_map.find(temp_q);
                        solver->addClause(~itr_p->second,~itr_q->second);
                    }
                }
                for (int i = counter_i; i <n; i ++){
                    int_pairs temp_p = make_pair(counter_i, counter_j);
                    itr_p = my_map.find(temp_p);
                    for (int j = 1; j < counter_j; j++){
                        int_pairs temp_q = make_pair(i, j);
                        itr_q = my_map.find(temp_q);
                        solver->addClause(~itr_p->second,~itr_q->second);
                    }
                }
            }
        }
        //clause 5 end

        bool res = solver->solve();
//        sleep(k*0.25);
        if (res == 1){
            for (int counter_i = 0; counter_i < n; counter_i++){
                for (int counter_k = 1; counter_k <= k; counter_k++){
                    int_pairs temp_j = make_pair(counter_i, counter_k);
                    itr_p = my_map.find(temp_j);
                    if(toInt(solver ->modelValue(itr_p->second)) == 0){
//                        cout<<counter_i<<" "<<counter_k<<endl;
                        result.push_back(counter_i);
                    }
                }
            }
            sort(result.begin(),result.end());
            break;
        }
        else{
            solver.reset(new Minisat::Solver());
            k++;
        }
    }
    return result;
}

int sort_vector(const Edge_list& temp_list){
    Edge_list my_list = temp_list;
    vector<int_pairs> vertex_vec;

    while (!my_list.empty()){
        bool start_found = false;
        bool end_found = false;
        Edge temp_edge = my_list.pop();//cout<<temp_edge.start<<", "<<temp_edge.end<<endl;
        for (auto & i : vertex_vec){
            if(i.first == temp_edge.start){
                i.second++;
                start_found = true;
            }
            else if (i.first == temp_edge.end){
                i.second ++;
                end_found = true;
            }
        }
        if(!start_found){
            vertex_vec.emplace_back(temp_edge.start,1);
        }
        if(!end_found){
            vertex_vec.emplace_back(temp_edge.end,1);
        }
    }
//    my_list.~Edge_list();

    sort(vertex_vec.begin(), vertex_vec.end(),
         [](const int_pairs &x, const int_pairs y) -> int {
             return x.second > y.second;
         });

    return vertex_vec[0].first;
}

vector<int> approx_vc_1(const Edge_list& temp_list) {

    Edge_list my_list = temp_list;
    vector<int> result;
    int high_degree_vertex;
    while (!my_list.empty()) {

        high_degree_vertex = sort_vector(my_list);
        for (int counter2 = 0; counter2 < my_list.size(); counter2++) {
            Edge temp_edge = my_list.pop();
            if (temp_edge.start == high_degree_vertex || temp_edge.end == high_degree_vertex) {
                result.push_back(high_degree_vertex);
                counter2--;
            } else {
                my_list.push(temp_edge);
            }
        }
    }

    sort(result.begin(), result.end());
    result.erase(unique(result.begin(), result.end()), result.end());//get the unique value from result;
    return result;
}

vector<int> approx_vc_2(const Edge_list& temp_list) {
    Edge_list my_list = temp_list;
    int my_length = my_list.get_v();
    vector<vector<int>> numEdge(my_length);
    vector<int> visited(my_length, 0);
    while (!my_list.empty()) {
        Edge temp_edge = my_list.pop();
        numEdge[temp_edge.start].push_back(temp_edge.end);
        numEdge[temp_edge.end].push_back(temp_edge.start);
    }
//    my_list.~Edge_list();

    for (int u = 0; u < my_length; ++u) {
        if (visited[u] == 0) {
            for (auto iter = numEdge[u].begin(); iter != numEdge[u].end(); ++iter) {
                if (visited[*iter] == 0) {
                    visited[*iter] = 1;
                    visited[u] = 1;
                    break;
                }
            }
        }
    }
    return visited;
}

void print_path(int start ,int end, vector<int>pre){
    if (start == end){
        cout<<start;
        return;
    }
    print_path(start,pre[end],pre);
    cout<<"-"<<end;
}

void draw_graph(const Edge_list& temp_list, int my_size, int start_point, int end_point)
{
    Edge_list my_list = temp_list;
    vector<bool> visit(my_size);
    vector<int> distance(my_size);
    vector<int> pre(my_size);
    try {
        if (start_point >= my_size || end_point >= my_size)
        {
            throw std::invalid_argument("Error: the vertex does not exist");
        }
        else if (start_point == end_point){
            cout<<start_point<<"-"<<end_point;
        }
        else{
            int my_graph[my_size][my_size];// maximum the graph
            for (int counter1 = 0; counter1 < my_size; counter1++){ //should be ok to take off the brackets here!!!!
                for(int counter2 = 0; counter2 < my_size; counter2++){
                    if (counter2 == counter1){
                        my_graph[counter1][counter2]= 0; // distance between the same vertex is 0;
                    }
                    else{
                        my_graph[counter1][counter2]= INF; // if the number is <0 then the path does not exist;
                    }
                }
            }
            int temp_size = my_list.size();
            for (int counter = 0; counter < temp_size; counter ++){ // change path by going through my_list
                Edge temp_edge = my_list.pop();
                my_graph[temp_edge.start][temp_edge.end] = 1;
                my_graph[temp_edge.end][temp_edge.start] = 1;
            }
            try {
                fill(distance.begin(), distance.end(), INF);
                fill(visit.begin(),visit.end(),false);
                for (int i = 0; i < my_size; ++i) {
                    pre[i] = i;
                }

                int u = -1; //find smallest u
                distance[start_point] = 0;
                for (int i = 0; i < my_size; ++i) {
                    int MIN = INF;                                  //keep distance
                    for (int j = 0; j < my_size; ++j)                     //start find next d[u]
                    {
                        if (!visit[j] && distance[j] < MIN) {
                            u = j;
                            MIN = distance[j];
                        }
                    }
                    if (u != -1)
                        visit[u] = true;                                 //mark visited
                    for (int v = 0; v < my_size; ++v) {
                        //go through list and find short path from point v
                        if (!visit[v] && distance[u] + my_graph[u][v] < distance[v]) {
                            distance[v] = distance[u] + my_graph[u][v];             //update shortest distance
                            pre[v] = u;                        //add u to path
                        }
                    }
                }

                if (!visit[end_point]){ // no path found
                    throw std::invalid_argument("Error: the path does not exist");
                }
                else{
                    print_path(start_point, end_point, pre);
                    cout<<endl;
                }

                for (int counter1 = 0; counter1 < my_size; counter1++) //should be ok to take off the brackets here!!!!
                    for (int counter2 = 0; counter2 < my_size; counter2++)
                        my_graph[counter1][counter2] = INF; // if the number is <0 then the path does not exist;
            }
            catch (const std::invalid_argument &e) {
                cerr << e.what()<<endl;
            }
        }
        my_list.~Edge_list();
        visit.clear();
        distance.clear();
        pre.clear();
    }
    catch (const std::invalid_argument &e) {
        cerr << e.what()<<endl;
    }
}

void print_result(thread_data my_data){

    cout<<"CNF-SAT-VC: ";
    if(my_data.cnf_result.empty()){
        cout<< "Exceed max wait time."<<endl;
    }
    else{
        for (unsigned long counter_1 = 0; counter_1<my_data.cnf_result.size()-1; counter_1++){
            cout<<my_data.cnf_result[counter_1]<<",";
        }
        cout<<my_data.cnf_result.back()<<endl;
//        cout<<"cnf running-time: "<<my_data.cnf_time<<" (us)."<<endl;
    }
    cout<< "APPROX-VC-1: ";
    for(unsigned long counter = 0; counter< my_data.approx_vc_1_result.size()-1;counter++)
    {
        cout << my_data.approx_vc_1_result[counter] << ",";
    }
    cout<<my_data.approx_vc_1_result.back()<<endl;
//    cout<<"Approx1 running-time: "<<my_data.vc_1_time<<" (us)."<<endl;

    cout<< "APPROX-VC-2: ";
    vector<int> temp_result2;
    for (int i = 0; i <= my_data.approx_vc_2_result.size(); ++i){
        if (my_data.approx_vc_2_result[i] == 1)
            temp_result2.push_back(i);
    }
    for (int i = 0; i < temp_result2.size()-1; i++){
        cout<<temp_result2[i]<<",";
    }
    cout << temp_result2.back()<<endl;
//    cout<<"Approx2 running-time: "<<my_data.vc_2_time<<" (us)."<<endl;
//    //reset my_data
    my_data.cnf_result.clear();
    my_data.approx_vc_1_result.clear();
    my_data.approx_vc_2_result.clear();
    my_data.cnf_time = 0;
    my_data.vc_1_time = 0;
    my_data.vc_2_time = 0;
//    my_data.thread_edge_list->~Edge_list();
}

// this thread controls I/O
void *thread_zero(void* ptr){
    // when thread 0 is on the run, lock thread 234
    sem_wait(&sem0);

    auto* new_ptr = (thread_data *) ptr;
    Edge_list new_list(1);//initialize a list
    string raw_input;
    regex pattern_V("[0-9]+"); //parse pattern for commend V
    regex pattern_E("[\\d]+"); //parse pattern for commend E
    smatch parsed_input_V;
    smatch parsed_input_E;
    vector <int> temp_edge_list; // this vertex holds all X,Y in input {<X1,Y1>....<Xn,Yn>}. X & Y may appear more than once.
    int my_length = 0;
    vector <int> temp_s_list;

//        pthread_mutex_lock(&work_mutex);
    while (getline(cin,raw_input)){
        if(raw_input.empty()){
            continue;
        }
        else if (strncmp(&raw_input.at(0), "V", 1) == 0) {
            try {
                regex_search(raw_input, parsed_input_V, pattern_V);
                if (parsed_input_V.length() == 0) {
                    throw std::invalid_argument("Error: V invalid entry format");
                } else {
                    auto sub = parsed_input_V.begin();// get the first index in parsed_input[]
                    stringstream temp_ss(sub->str()); // covert string to sstream
                    my_length = 0;
                    temp_ss >> my_length; // covert sstream to int
                    new_list.~Edge_list();
                    new_list = Edge_list(my_length);
                    new_ptr -> thread_edge_list = &new_list;
                }
            }
            catch (const std::invalid_argument &e) {
                cerr << e.what() << endl;
            }
        }
        else if (strncmp(&raw_input.at(0), "E", 1) == 0) {
            try {
                string::const_iterator iterStart = raw_input.begin();
                string::const_iterator iterEnd = raw_input.end();
                string temp;
                while (regex_search(iterStart, iterEnd, parsed_input_E, pattern_E)) {
                    temp = parsed_input_E[0];
                    stringstream temp_ss2(temp);
                    int temp_number = 0;
                    temp_ss2 >> temp_number; // covert sstream to int
                    if (temp_number >= new_ptr->thread_edge_list->get_v()){
                        // when this error is throw, clean the vector.
                        throw std::invalid_argument("Error: If the integer that follows the V is i then we assume that the vertices are identified by 0,...,i-1");
                    }
                    temp_edge_list.push_back(temp_number); // covert sstream to int
                    iterStart = parsed_input_E[0].second;    //更新搜索起始位置,搜索剩下的字符串
                }
                for (unsigned long counter = 0; counter<temp_edge_list.size();counter = counter+2){
                    Edge my_edge{};
                    my_edge.start = temp_edge_list[counter];
                    my_edge.end = temp_edge_list[counter+1];
                    new_list.push(my_edge);
                }
            }
            catch (const std::invalid_argument& e) {
                cerr << e.what() << endl;
                temp_edge_list.clear();
            }

            temp_edge_list.clear();
        }
        else if (strncmp(&raw_input.at(0), "s", 1) == 0){
            stringstream ss_s(raw_input);
            /* Running loop till the end of the stream */
            string temp;
            int check_int = 0;
            while (!ss_s.eof()) {
                /* extracting word by word from stream */
                ss_s >> temp;
                /* Checking the given word is integer or not */
                if (stringstream(temp) >> check_int)
                    temp_s_list.push_back(check_int);
                /* To save from space at the end of string */
                temp = "";
            }
            try{
                draw_graph(*(new_ptr->thread_edge_list),new_ptr->thread_edge_list->get_v(),temp_s_list[0],temp_s_list[1]);
                temp_s_list.clear();
            }
            catch (const std::invalid_argument& e) {
                cerr<< e.what() << endl;
                temp_s_list.clear();
            }
        }
        else {
//            new_ptr->thread_edge_list->~Edge_list();
            continue;
        }

        // unlock thread 1,2,3
        if(!new_ptr->thread_edge_list->empty()){
            sem_post(&sem1);
            sem_post(&sem2);
            sem_post(&sem3);

            sem_wait(&sem4);
            //ready to print;

            p_mutex.lock();
            print_result(*new_ptr);
            p_mutex.unlock();
        }
        if(new_ptr->thread_edge_list == NULL){
            //unlock input
            sem_post(&sem0);
        }
    }
    return nullptr;
}

void *thread_one(void *ptr){
    sem_wait(&sem1); //lock sem1
    auto* new_ptr1 = (thread_data *) ptr;
    Edge_list temp(new_ptr1->thread_edge_list->size());
    temp = *(new_ptr1->thread_edge_list);
    if(!temp.empty()){
        start1 = high_resolution_clock::now();
        // pass in the thread list's value into getCNF(), then store its result int thread data
        new_ptr1->cnf_result = getCNF(temp);
        end1 = high_resolution_clock::now();
        new_ptr1->cnf_time = duration_cast<microseconds>(end1 - start1).count();
    }

    return nullptr;
}

void *thread_two(void *ptr){
    sem_wait(&sem2); //lock sem2
    auto* new_ptr2 = (thread_data *) ptr;

    Edge_list temp(new_ptr2->thread_edge_list->size());
    temp = *(new_ptr2->thread_edge_list);
    if(!temp.empty()){
        start1 = high_resolution_clock::now();
        new_ptr2->approx_vc_1_result = approx_vc_1(temp);
        end1 = high_resolution_clock::now();
        new_ptr2->vc_1_time = duration_cast<microseconds>(end1 - start1).count();
    }
//    sem_post(&sem1);
    return nullptr;
}

void *thread_three(void *ptr){
    sem_wait(&sem3); //lock se3
    auto* new_ptr3 = (thread_data *) ptr;
    Edge_list temp(new_ptr3->thread_edge_list->size());
    temp = *(new_ptr3->thread_edge_list);
    if(!(temp.empty())){
        start1 = high_resolution_clock::now();
        new_ptr3->approx_vc_2_result = approx_vc_2(temp);
        end1 = high_resolution_clock::now();
        new_ptr3->vc_2_time = duration_cast<microseconds>(end1 - start1).count();
    }
//    sem_post(&sem2);

    return nullptr;
}

int main() {

    thread_data my_data;//first initialize

    int ret;

    pthread_t thread_0;// io
    pthread_t thread_1;// cnf
    pthread_t thread_2;// method 1
    pthread_t thread_3;// method 2

    int sem[5];
    sem[0] = sem_init(&sem0,0,1); // inti sem1 = 1; control input
    sem[1] = sem_init(&sem1,0,0); // inti sem2 = 0;
    sem[2] = sem_init(&sem2,0,0); // inti sem3 = 0;
    sem[3] = sem_init(&sem3,0,0); // inti sem4 = 0;
    sem[4] = sem_init(&sem4,0,0); // control output

    for (int counter = 0; counter < 4; counter ++ ){
        if(sem[counter] != 0){
            cerr<<"Error: failed to create sem "<<counter<<"."<<endl;
        }
    }

    ret = pthread_create(&thread_0,nullptr,thread_zero, (void *) &my_data);//dereference input_list because thread_one take in a pointer as argument
    if(ret != 0){
        cerr<<"Error: failed to create thread 1"<<endl;
    }

    while (true){
        ret = pthread_create(&thread_1,nullptr,thread_one, (void *)&my_data);//dereference input_list because thread_one take in a pointer as argument
        if(ret != 0){
            cerr<<"Error: failed to create thread 2"<<endl;
            break;
        }
        ret = pthread_create(&thread_2,nullptr,thread_two, (void *)&my_data);//dereference input_list because thread_one take in a pointer as argument
        if(ret != 0){
            cerr<<"Error: failed to create thread 3"<<endl;
            break;
        }
        ret = pthread_create(&thread_3,nullptr,thread_three, (void *)&my_data);//dereference input_list because thread_one take in a pointer as argument
        if(ret != 0){
            cerr<<"Error: failed to create thread 4"<<endl;
            break;
        }


        pthread_join(thread_1, nullptr);
        pthread_join(thread_2, nullptr);
        pthread_join(thread_3, nullptr);
        sem_post(&sem4);//unlock output;
    }


    pthread_join(thread_0, nullptr);
    return 0;
}
