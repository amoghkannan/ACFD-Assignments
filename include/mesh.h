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
        void operator*=(double val){ x*=val;y*=val;};
        Node operator/(double val){ return Node(x/val,y/val);};
        void operator/=(double val){ x/=val;y/=val;};
        Node operator-(Node& n){return Node(x-n.x,y-n.y);};
        Node operator+(Node& n){return Node(x+n.x,y+n.y);};
        wp norm2(){return sqrt(pow(x,2.0)+pow(y,2.0));};
};

struct Vec2{
        Node st,en;

        Vec2(){};
        Vec2(Node nd1,Node nd2): st(nd1), en(nd2) {};
        Vec2(wp xComp,wp yComp): st(Node(0.0,0.0)), en(Node(xComp,yComp)) {};
        Vec2 operator-(){return Vec2(en,st);};
        Vec2 operator-(Vec2 otherVec){return Vec2(st-otherVec.st,en-otherVec.en);};
        void operator*=(double val){st*=val;en*=val;};
        Vec2 operator*(double val){return Vec2(st*val,en*val);};
        void operator/=(double val){st/=val;en/=val;};
        Vec2 operator/(double val){return Vec2(st/val,en/val);};
        wp dotProduct(Vec2 otherVec){return (en.x-st.x)*(otherVec.en.x-otherVec.st.x)+
                                            (en.y-st.y)*(otherVec.en.y-otherVec.st.y);};
        wp norm2(){return (en-st).norm2();};
};

class Mesh: public Grid<Node>{

private:

public:

Grid<wp>volumes;
Grid<Vec2>normalsI;
Grid<Vec2>normalsJ;
static wp calcVolume(std::vector<Node>nodes);
void calcVolumes();
void calcAreas();
Mesh(const Mesh& otherMesh): Grid<Node>(otherMesh){};
Mesh(Box box_, std::array<idtype,2>indexType_): Grid<Node>(box_, indexType_) {};
void setGhostNodes();
Mesh coarsen();
void operator+(Mesh& otherMesh)=delete;
void operator*(Mesh& otherMesh)=delete;
void operator/(Mesh& otherMesh)=delete;
Node cc(int i, int j);
static wp pDist(Vec2 line, Node point, Vec2 normal);
~Mesh();

};
