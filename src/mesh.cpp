#include"mesh.h"

Mesh::Mesh(int imx,int jmx,int bufW, int bufE, int bufS, int bufN){
        this->imx=imx;
        this->jmx=jmx;
        this->bufE=bufE;
        this->bufW=bufW;
        this->bufN=bufN;
        this->bufS=bufS;

        data=new Node[(imx+bufE+bufW)*(jmx+bufN+bufS)];

        logger.log("Debug:Mesh initialization",1);
};

std::array<int,6> Mesh::size(){
        std::array<int,6>ans={imx,jmx,bufW,bufE,bufS,bufN};
        return ans;
};

Node& Mesh::operator()(int i,int j){
        
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

Mesh::~Mesh(){
delete[] data;
logger.log("Debug: destroying mesh",1);
};
