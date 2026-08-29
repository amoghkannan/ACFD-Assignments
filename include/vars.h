#pragma once
#include"utils.h"
#include<array>

class Grid{

private:

int imx=0;
int jmx=0;
int bufE,bufW,bufN,bufS=0;

wp* data=nullptr;

public:

Grid(const Grid& otherGrid){
       
       imx=otherGrid.imx;
       jmx=otherGrid.jmx;
       bufE=otherGrid.bufE;
       bufW=otherGrid.bufW;
       bufN=otherGrid.bufN;
       bufS=otherGrid.bufS;

       data=new wp[(imx+bufE+bufW)*(jmx+bufN+bufS)];
        
       for(int j=0;j<jmx+bufN+bufS;j++){
                for(int i=0;i<imx+bufE+bufW;i++){
                        data[j*(imx+bufE+bufW)+i]=otherGrid.data[j*(imx+bufE+bufW)+i];
                };
        };

};

Grid(int imx, int jmx, int bufW, int bufE, int bufS, int bufN){
        this->imx=imx;
        this->jmx=jmx;
        this->bufE=bufE;
        this->bufW=bufW;
        this->bufN=bufN;
        this->bufS=bufS;

        data=new wp[(imx+bufE+bufW)*(jmx+bufN+bufS)];

};

std::array<int,6>size();

wp& operator()(int i, int j);
void operator=(Grid& otherGrid);
void operator+(Grid& otherGrid);
void operator*(wp val);
void operator/(wp val);

~Grid();

};
