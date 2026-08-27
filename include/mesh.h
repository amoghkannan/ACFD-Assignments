#pragma once
#include"utils.h"
#include<stdlib.h>
#include<array>

struct Node{
        wp x;
        wp y;
};

class Mesh{

private:

int imx=0;
int jmx=0;
int bufE,bufW,bufN,bufS=0;

Node* data;

public:

Mesh(){
};

Mesh(int imx, int jmx, int bufW, int bufE, int bufS, int bufN);

Mesh(const Mesh& meshIn){
        Mesh(meshIn.imx,meshIn.jmx,meshIn.bufW,meshIn.bufE,meshIn.bufS,meshIn.bufN);
        *(this->data)=*(meshIn.data);
};

std::array<int,6>size();

Node& operator()(int i, int j);

~Mesh();

};
