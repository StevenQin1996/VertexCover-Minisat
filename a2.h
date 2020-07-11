//
// Created by Shiyun Qin on 2019-09-30.
//

#ifndef ECE650_A2_A2_H
#define ECE650_A2_A2_H
//
// Created by Shiyun Qin on 2019-09-30.
//
#ifndef ECE650_A2_ECE650A2_H
#define ECE650_A2_ECE650A2_H

using namespace std;

struct Edge{
    int start;
    int end;
};

class Edge_list{
    //need to implement a way to output Edge_list
    friend std::ostream& operator<<(std::ostream&, const Edge_list&);

private:
    int length;
    int head;
    int v;//size
    int tail;
    Edge *data;// data type should be Edge
//    vector<int> vector_temp_edge_list;

public:
    explicit Edge_list(int i); //constructor //study explicit
    virtual ~Edge_list(); //destroyer
    int size()const; //get size
    bool empty()const;
    int get_v() const;
//    void set_unique_vertex_list(vector<int> blee_list) {vector_temp_edge_list = blee_list;}
//    vector<int> get_unique_vertex_list() const;
    void push(Edge); // add a edge
    Edge pop();
};

//initialize Edge_list with the size of user input
Edge_list::Edge_list(int i) {
    length = i * i;
    data = new Edge[length];
    head = 0;
    tail = 0;
    v = i;
}

//destoryer
Edge_list::~Edge_list() {
    data = nullptr;
    delete[] data;
}

//get the size
int Edge_list::size() const {
    int count=0;
    if(tail>=head)
        count = tail-head;
    else
        count = tail+ (length - head);
    return count;
}

//check whether the Edge_list is empty
bool Edge_list::empty()const
{
    return (size()==0);
}

//vector<int> Edge_list::get_unique_vertex_list() const{
//    vector <int> unique_vertex_list; //this vertex is sorted from temp_edge_list。use set temp_set
//    set<int> temp_set(vector_temp_edge_list.begin(), vector_temp_edge_list.end());
//    unique_vertex_list.assign(temp_set.begin(),temp_set.end());
//    return unique_vertex_list;
//}

//push Edge into Edge_list
void Edge_list::push(Edge my_Edge) {
    try {
        if (my_Edge.start != my_Edge.end){
            if (size() == length-1) // if out of size, then create a larger list
            {
                Edge* larger_list = new Edge[length*2];

                for (int x=0; x < length; x++)
                {
                    larger_list[x]= data[head++];

                    if (head == length)
                        head=0;
                }
                tail = length-1;
                head = 0;
                length = length * 2;

                delete [] data;
                data = larger_list;
            }

            data[tail]= my_Edge;
            tail++;

            if (tail==length )
            {
                tail=0;
            }

//            cout<< "push done!!!\n";
        }
        else{
            throw std::invalid_argument("Edge can not have the same start and end\n");
        }
    }
    catch (const char* e){
        cerr<<"Error:"<<e;
    }
}

Edge Edge_list::pop()//pop
{
    Edge first_out=data[head];
    head++;
    if (head==length)
    {
        head=0;
    }
    return first_out;
}

std::ostream& operator<<(std::ostream& lhs , const Edge_list& rhs)//name ostream and print out IntegerQueue.
{
    lhs<<"{";
    if(!rhs.empty())//if rhs is empty, there should not be any number in side the {};
    {
        int counter = 0;
        int head = rhs.head;
        while (counter < rhs.size()-1){

            if(head ==rhs.length)
                head = 0;
            lhs << "<"<<rhs.data[head].start<<","<<rhs.data[head].end<<">, ";
            head++;
            counter ++;
        }
        lhs << "<"<<rhs.data[head].start<<","<<rhs.data[head].end<<">";
    }
    lhs<<"}";
    return lhs;
}

int Edge_list::get_v() const{
    return v;
}


#endif //ECE650_A2_ECE650A2_H

#endif //ECE650_A2_A2_H
