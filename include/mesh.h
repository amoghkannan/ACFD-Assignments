#pragma once
#include"utils.h"
#include<stdlib.h>
#include<array>
#include<vars.h>

struct Node{
        wp x;
        wp y;

        Node(): x(0.0),y(0.0){};
        Node(wp x_, wp y_): x(x_), y(y_){};
        Node operator*(double val){ return Node(x*val,y*val);};
        Node operator-(Node& n){return Node(x-n.x,y-n.y);};
        Node operator+(Node& n){return Node(x+n.x,y+n.y);};
        wp norm2(){return sqrt(pow(x,2.0)+pow(y,2.0));};
};

struct Vec2{
        Node st,en;

        Vec2(){};
        Vec2(Node nd1,Node nd2): st(nd1), en(nd2) {};
        wp norm2(){return (en-st).norm2();};
};

class Mesh: public Grid<Node>{

private:

public:

Grid<wp>volumes;
static wp calcVolume(std::vector<Node>nodes);
void calcVolumes();
Mesh(const Mesh& otherMesh): Grid<Node>(otherMesh){};
Mesh(Box box_, std::array<idtype,2>indexType_): Grid<Node>(box_, indexType_) {};
void setGhostNodes();
Mesh coarsen();
void operator+(Mesh& otherMesh)=delete;
void operator*(Mesh& otherMesh)=delete;
void operator/(Mesh& otherMesh)=delete;
~Mesh();

};
