#include"vars.h"

std::array<int,6> Grid::size(){
        std::array<int,6>ans={imx,jmx,bufW,bufE,bufS,bufN};
        return ans;
};

wp& Grid::operator()(int i,int j){
        
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

Grid::~Grid(){

        delete[] data;
};
