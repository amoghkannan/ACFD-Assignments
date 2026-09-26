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
Mesh(Box box_, std::array<idtype,2>indexType_): Grid<Node>(box_, indexType_) {};

void operator+(Mesh& otherMesh)=delete;
void operator*(Mesh& otherMesh)=delete;
void operator/(Mesh& otherMesh)=delete;

~Mesh();

};
