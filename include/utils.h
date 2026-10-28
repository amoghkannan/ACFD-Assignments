#pragma once
#include"logger.h"
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<utility>
#include<iomanip>
#include<chrono>
#include<algorithm>
#include<unordered_map>
#include<tuple>

#define wp double
#define BIG 1E9
#define SMALL 1E-14

using Point=std::pair<int,int> ;
using boundMatEntry=std::pair<Point,wp>;
using boundMatRow=std::vector<boundMatEntry>;

extern Logger logger;

#ifdef Debug
#define LOGCODE 1
#else
#define LOGCODE 2
#endif

wp sign(wp in);

template<typename key, typename val>
class dictionary{

private:

std::vector<key> keys;

std::unordered_map<key,val>data;

public:

val& operator[](key keyIn){

        auto it=data.find(keyIn);
        if(it==data.end()){
                keys.push_back(keyIn);
        };

        return data[keyIn];
};

val& operator()(int ind){

        return data[keys[ind]]; 

};

bool hasKey(key keyIn){
       
        auto it=data.find(keyIn);
        if(it==data.end()){
                return false;
        };

        return true;
};

std::vector<key>getKeys(){
        return keys;
};

};
