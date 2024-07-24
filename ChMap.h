#ifndef CH_MAP_H_093002172024
#define CH_MAP_H_093002172024

#include <map>
#include <string>
using namespace std;

class ChMap{
public:
    // constructor and destructor
    ChMap();
    ~ChMap();

public:
    void Set(std::vector<string>  strips, int n_of_ch); // set a ch map. The argument is like "{"8R", "", "2L", "9R", ""}".
    int Get(string stripSide) const; // return a discri ch corresponding to a given stripside(e.g. "8R") in the ch map

    //void Set(string strip, int ch); オーバーライド？する
/*
private:
    void CheckKey(string stripside[]);
*/

private:
    std::map<string, int> m_chMap;


};


ChMap::ChMap(){
    std::map<string, int> m_chMap{};
}

ChMap::~ChMap(){
    // メンバーにポインタがいないので何もしなくてよし(ホンマか？)
}

void ChMap::Set(std::vector<string> strips, int n_of_ch){
    for(int i=0; i < n_of_ch; i++){
        m_chMap[strips[i]] = i;
    }
}

int ChMap::Get(string stripSide) const {
    return m_chMap.at(stripSide);
}
#endif