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

Grid(){};

Grid<T>(const Grid<T>&otherGrid){
       
       imx=otherGrid.imx;
       jmx=otherGrid.jmx;
       bufE=otherGrid.bufE;
       bufW=otherGrid.bufW;
       bufN=otherGrid.bufN;
       bufS=otherGrid.bufS;

       if(data!=nullptr) delete [] data;
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

        if(data!=nullptr) delete [] data;
        data=new T[(imx+bufE+bufW)*(jmx+bufN+bufS)];

};

std::array<int,6>size();

T& operator()(int i, int j);
void operator=(Grid otherGrid);
void operator+=(Grid& otherGrid);
void operator-=(Grid& otherGrid);
void operator-();
void operator*=(wp val);
Grid<T> operator*(wp val);
void operator/=(wp val);
void print(std::ostream& os);
wp norm2();
wp dotProduct(Grid& otherGrid);
void initVal(T val);
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
void Grid<T>::operator=(Grid<T> otherGrid){

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

                if(data!=nullptr) delete [] data;
                data=new T[(imx+bufE+bufW)*(jmx+bufN+bufS)];
          };


        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] = otherGrid(i,j);
                };
        };
};

template<typename T>
void Grid<T>::operator+=(Grid<T>& otherGrid){
        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] = data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] + otherGrid(i,j);
                };
        };
};

template<typename T>
void Grid<T>::operator-=(Grid<T>& otherGrid){
        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] = data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] - otherGrid(i,j);
                };
        };
};

template<typename T>
void Grid<T>::operator-(){
        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] = -data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1];
                };
        };
};

template<typename T>
void Grid<T>::operator*=(wp val){
        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] = data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] * val;
                };
        };
};

template<typename T>
Grid<T> Grid<T>::operator*(wp val){

        Grid<T> newGrid(*this);

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        newGrid.data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1]=
                        data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] * val;
                };
        };
};

template<typename T>
void Grid<T>::operator/=(wp val){
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
wp Grid<T>::norm2(){
        wp ans=0.0;

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        ans=ans+pow(data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1],2.0);
                };
        };

        ans=sqrt(ans);
        return ans;
};

template<typename T>
wp Grid<T>::dotProduct(Grid<T>& otherGrid){
        wp ans=0.0;

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        ans=ans+data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1]*
                      otherGrid.data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1];
                };
        };

        return ans;
};

template<typename T>
void Grid<T>::initVal(T val){
        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1]=val;
                };
        };

};

template<typename T>
Grid<T>::~Grid(){

        delete[] data;

};
