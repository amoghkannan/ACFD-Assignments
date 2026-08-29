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

void Grid::operator=(Grid& otherGrid){
        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] = otherGrid(i,j);
                };
        };
};

void Grid::operator+(Grid& otherGrid){
        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] = data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] + otherGrid(i,j);
                };
        };
};

void Grid::operator*(wp val){
        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] = data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] * val;
                };
        };
};

void Grid::operator/(wp val){
        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] = data[(j+bufS-1)*(imx+bufE+bufW)+i+bufW-1] / val;
                };
        };
};

Grid::~Grid(){

        delete[] data;

        logger.log("Debug: Destroying grid",1);
};
