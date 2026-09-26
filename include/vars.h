#pragma once
#include"utils.h"
#include<array>

enum idtype{
        NODE,
        CELL
};

struct Box{ //(Nodal) box in index space on which the grid is defined
        int imx=0; //Based on nodes; not representative of actual number of I items
        int jmx=0; //Ditto
        int bufW,bufE,bufS,bufN=0;
        Box(){};
        Box(int imx_, int jmx_, int bufW_, int bufE_, int bufS_, int bufN_): imx(imx_),jmx(jmx_),bufW(bufW_),
        bufE(bufE_),bufS(bufS_),bufN(bufN_){};
        bool operator!=(Box& otherBox){
                return (imx!=otherBox.imx) || (jmx!=otherBox.jmx) || (bufE!=otherBox.bufE) || (bufW!=otherBox.bufW) 
                       || (bufN!=otherBox.bufN) || (bufS!=otherBox.bufS);
        };
};

template<typename T>
class Grid{

protected:

T* data=nullptr;

public:

Box box;
std::array<idtype,2> indexType; //Node, Iface, Jface, cell
std::array<int,2> intVect; //Actual dimension: eg 1:imx-1 for cells
int trueSize=0; //Grid sometimes allocated more space than it needs to save on allocation time in future
//Even then, indexing into grid done using intVect, buffers 

//Copying rules:
//Hidden cells are never copied
//If current storage is larger than copied grid, current storage is retained
//Hidden data never copied

Grid(){};

Grid<T>(const Grid<T>&otherGrid){
      
       if(data!=nullptr){
                if((otherGrid.intVect[0]+otherGrid.box.bufE+otherGrid.box.bufW)*
                   (otherGrid.intVect[1]+otherGrid.box.bufN+otherGrid.box.bufS)>trueSize){ 
                        delete [] data;
                        trueSize=(otherGrid.intVect[0]+otherGrid.box.bufE+otherGrid.box.bufW)*
                           (otherGrid.intVect[1]+otherGrid.box.bufN+otherGrid.box.bufS);
                        data=new T[trueSize];
                };

       }
       else{
                trueSize=(otherGrid.intVect[0]+otherGrid.box.bufE+otherGrid.box.bufW)*
                           (otherGrid.intVect[1]+otherGrid.box.bufN+otherGrid.box.bufS);
                data=new T[trueSize];
       };

       box=otherGrid.box;

       indexType=otherGrid.indexType;
       intVect=otherGrid.intVect;

       for(int j=1;j<=intVect[1];j++){
                for(int i=1;i<=intVect[0];i++){
                        (*this)(i,j)=otherGrid.data[(j+box.bufS-1)*(intVect[0]+box.bufE+box.bufW)+i+box.bufW-1];
                };
        };

};

Grid(Box box_, std::array<idtype,2>indexType_){

        box=box_;
        indexType=indexType_;
        intVect[0]=indexType[0]==CELL?box.imx-1:box.imx;
        intVect[1]=indexType[1]==CELL?box.jmx-1:box.jmx;

       if(data!=nullptr){
                if((intVect[0]+box_.bufE+box_.bufW)*(intVect[1]+box_.bufN+box_.bufS)>trueSize){ 
                        delete [] data;
                        trueSize=(intVect[0]+box_.bufE+box_.bufW)*(intVect[1]+box_.bufN+box_.bufS);
                        data=new T[trueSize];
                };

       }
       else{
               trueSize=(intVect[0]+box_.bufE+box_.bufW)*(intVect[1]+box_.bufN+box_.bufS);
               data=new T[trueSize];
       };

};

std::array<int,6>size();

T& operator()(int i, int j);
void operator=(const Grid& otherGrid);
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
        std::array<int,6>ans={intVect[0],intVect[1],box.bufW,box.bufE,box.bufS,box.bufN};
        return ans;
};

template<typename T>
T& Grid<T>::operator()(int i,int j){
        
        if(i<1-box.bufW || i>intVect[0]+box.bufE){
                std::cout<<"Invalid I index, exiting";
                std::exit(-1);
        };

        if(j<1-box.bufS || j>intVect[1]+box.bufN){
                std::cout<<"Invalid J index, exiting";
                std::exit(-1);
        };

        return data[(j+box.bufS-1)*(intVect[0]+box.bufE+box.bufW)+i+box.bufW-1];
};

template<typename T>
void Grid<T>::operator=(const Grid<T>& otherGrid){
       if(data!=nullptr){
                if((otherGrid.intVect[0]+otherGrid.box.bufE+otherGrid.box.bufW)*
                   (otherGrid.intVect[1]+otherGrid.box.bufN+otherGrid.box.bufS)>trueSize){ 
                        delete [] data;
                        trueSize=(otherGrid.intVect[0]+otherGrid.box.bufE+otherGrid.box.bufW)*
                           (otherGrid.intVect[1]+otherGrid.box.bufN+otherGrid.box.bufS);
                        data=new T[trueSize];
                };

       }
       else{
                trueSize=(otherGrid.intVect[0]+otherGrid.box.bufE+otherGrid.box.bufW)*
                           (otherGrid.intVect[1]+otherGrid.box.bufN+otherGrid.box.bufS);
                data=new T[trueSize];
       };

       box=otherGrid.box;

       indexType=otherGrid.indexType;
       intVect=otherGrid.intVect;

       for(int j=1;j<=intVect[1];j++){
                for(int i=1;i<=intVect[0];i++){
                        (*this)(i,j)=otherGrid.data[(j+box.bufS-1)*(intVect[0]+box.bufE+box.bufW)+i+box.bufW-1];
                };
        };

};

template<typename T>
void Grid<T>::operator+=(Grid<T>& otherGrid){
        if(box!=otherGrid.box || intVect!=otherGrid.intVect || indexType!=otherGrid.indexType)
        logger.log("Warning, adding non-equivalent fields!",1);

        for(int j=1;j<=intVect[1];j++){
                for(int i=1;i<=intVect[0];i++){
                        (*this)(i,j) = (*this)(i,j) + otherGrid(i,j);
                };
        };
};

template<typename T>
void Grid<T>::operator-=(Grid<T>& otherGrid){
        if(box!=otherGrid.box || intVect!=otherGrid.intVect || indexType!=otherGrid.indexType)
        logger.log("Warning, subtracting non-equivalent fields!",1);

        for(int j=1;j<=intVect[1];j++){
                for(int i=1;i<=intVect[0];i++){
                        (*this)(i,j) = (*this)(i,j) - otherGrid(i,j);
                };
        };
};

template<typename T>
void Grid<T>::operator-(){

        for(int j=1;j<=intVect[1];j++){
                for(int i=1;i<=intVect[0];i++){
                        (*this)(i,j) = -(*this)(i,j);
                };
        };
};

template<typename T>
void Grid<T>::operator*=(wp val){
        for(int j=1;j<=intVect[1];j++){
                for(int i=1;i<=intVect[0];i++){
                        (*this)(i,j) = (*this)(i,j) * val;
                };
        };
};

template<typename T>
Grid<T> Grid<T>::operator*(wp val){

        Grid<T> newGrid(*this);
        newGrid*=val;
        return newGrid;

};

template<typename T>
void Grid<T>::operator/=(wp val){
        if(val==0.0) logger.log("Error, grid division by 0",1);

        for(int j=1;j<=intVect[1];j++){
                for(int i=1;i<=intVect[0];i++){
                        (*this)(i,j) = (*this)(i,j) / val;
                };
        };
};

template<typename T>
void Grid<T>::print(std::ostream& os){
        for(int j=1;j<=intVect[1];j++){
                for(int i=1;i<=intVect[0];i++){
                        os<<i<<"\t"<<j<<"\t"<<(*this)(i,j)<<std::endl;
                };
        };
};

template<typename T>
wp Grid<T>::norm2(){
        wp ans=0.0;

        for(int j=1;j<=intVect[1];j++){
                for(int i=1;i<=intVect[0];i++){
                        ans=ans+pow((*this)(i,j),2.0);
                };
        };

        ans=sqrt(ans);
        return ans;
};

template<typename T>
wp Grid<T>::dotProduct(Grid<T>& otherGrid){
        if(box!=otherGrid.box || intVect!=otherGrid.intVect || indexType!=otherGrid.indexType)
        logger.log("Warning, dot product of non-equivalent fields!",1);

        wp ans=0.0;

        for(int j=1;j<=intVect[1];j++){
                for(int i=1;i<=intVect[0];i++){
                        ans=ans+(*this)(i,j)*otherGrid(i,j);
                };
        };

        return ans;
};

template<typename T>
void Grid<T>::initVal(T val){
        for(int j=1;j<=intVect[1];j++){
                for(int i=1;i<=intVect[0];i++){
                        (*this)(i,j)=val;
                };
        };

};

template<typename T>
Grid<T>::~Grid(){

        delete[] data;

};
