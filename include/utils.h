#pragma once
#include"logger.h"
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<utility>
#include<iomanip>

#define wp long double

using Point=std::pair<int,int> ;
using boundMatEntry=std::pair<Point,wp>;
using boundMatRow=std::vector<boundMatEntry>;

extern Logger logger;

#ifdef Debug
#define LOGCODE 1
#else
#define LOGCODE 2
#endif

template<typename key, typename val>
class dictionary{

private:

std::vector<key> keys;
std::vector<val> vals;

public:

val& operator[](key keyIn){

        int counter=0;
        int n=vals.size();
        
        for(int i=0; i<n;i++){
                key currKey=keys[i];

                if(currKey==keyIn){
                        return vals[counter];
                };

                counter=counter+1;
        };
        
        keys.resize(n+1);
        keys[n]=keyIn;
        vals.resize(n+1);
        return vals[n];

};


};
