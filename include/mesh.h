#pragma once
#include"utils.h"
#include<stdlib.h>
#include<array>
#include<vars.h>

struct Node{
        wp x;
        wp y;
};

class Mesh: public Grid<Node>{

private:

public:

Mesh(const Mesh& otherMesh): Grid<Node>(otherMesh){};
Mesh(int imx, int jmx, int bufW, int bufE, int bufS, int bufN): Grid<Node>(imx,jmx,bufW,bufE,bufS,bufN) {};

void operator+(Mesh& otherMesh)=delete;
void operator*(Mesh& otherMesh)=delete;
void operator/(Mesh& otherMesh)=delete;

~Mesh();

};
