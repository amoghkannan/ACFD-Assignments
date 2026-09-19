#pragma once
#include"utils.h"
#include<array>

template<typename T>
class Grid{

protected:

int imx=0;
int jmx=0;
int bufE,bufW,bufN,bufS=0;

T* data=nullptr;

public:

Grid(const Grid& otherGrid){
       
       imx=otherGrid.imx;
       jmx=otherGrid.jmx;
       bufE=otherGrid.bufE;
       bufW=otherGrid.bufW;
       bufN=otherGrid.bufN;
       bufS=otherGrid.bufS;

       data=new T[(imx+bufE+bufW)*(jmx+bufN+bufS)];
        
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

        data=new T[(imx+bufE+bufW)*(jmx+bufN+bufS)];

};

std::array<int,6>size();

T& operator()(int i, int j);
void operator=(Grid& otherGrid);
void operator+(Grid& otherGrid);
void operator*(wp val);
void operator/(wp val);
void print(std::ostream& os);
~Grid();

};

template<typename T>
std::array<int,6> Grid<T>::size(){
        std::array<int,6>ans={imx,jmx,bufW,bufE,bufS,bufN};
        return ans;
};

template<typename T>
T& Grid<T>::operator()(int i,int j){
        
        if(i<1-bufW || i>imx+bufE){
                std::cout<<"Invalid I index, exiting";
                std::exit(-1);
        };

        if(j<1-bufS || j>jmx+bufE){
                std::cout<<"Invalid J index, exiting";
                std::exit(-1);
        };

        return data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1];
};

template<typename T>
void Grid<T>::operator=(Grid<T>& otherGrid){

       if(imx!=otherGrid.imx ||
          jmx!=otherGrid.jmx ||
          bufE!=otherGrid.bufE ||
          bufW!=otherGrid.bufW ||
          bufN!=otherGrid.bufN ||
          bufS!=otherGrid.bufS){
                delete[] data;
                imx=otherGrid.imx;
                jmx=otherGrid.jmx;
                bufE=otherGrid.bufE;
                bufW=otherGrid.bufW;
                bufN=otherGrid.bufN;
                bufS=otherGrid.bufS;

                data=new T[(imx+bufE+bufW)*(jmx+bufN+bufS)];
          };


        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] = otherGrid(i,j);
                };
        };
};

template<typename T>
void Grid<T>::operator+(Grid<T>& otherGrid){
        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] = data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] + otherGrid(i,j);
                };
        };
};

template<typename T>
void Grid<T>::operator*(wp val){
        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] = data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] * val;
                };
        };
};

template<typename T>
void Grid<T>::operator/(wp val){
        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] = data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] / val;
                };
        };
};

template<typename T>
void Grid<T>::print(std::ostream& os){
        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        os<<i<<"\t"<<j<<"\t"<<data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1]<<std::endl;
                };
        };
};

template<typename T>
Grid<T>::~Grid(){

        delete[] data;

        logger.log("Debug: Destroying grid",1);
};
