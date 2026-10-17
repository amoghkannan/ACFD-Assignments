#include"utils.h"
#include"vars.h"
#include"mesh.h"

#define BIG 1E9
#define SMALL 1E-15

struct interfaceLine{
        Vec2 l;
        Vec2 n;
        interfaceLine(){};
        interfaceLine(Vec2 l_, Vec2 n_): l(l_),n(n_){};
};

struct IBDistStat{
        wp alpha;
        wp beta;
        wp dist;
        Node touchPoint;
        IBDistStat(){};
        IBDistStat(wp alpha_, wp beta_, wp dist_, Node touchPoint_):
        alpha(alpha_), beta(beta_), dist(dist_), touchPoint(touchPoint_) {};
        void operator=(IBDistStat other){
                alpha=other.alpha;
                beta=other.beta;
                dist=other.dist;
                touchPoint=other.touchPoint;
        };
};

class Interface{

private:
int nIB;

std::vector<interfaceLine>iblines;
Grid<wp>f; //Volume fraction function
Grid<wp>phi; //Level set function
Grid<wp>cellClass;
Mesh *mesh;

std::string meshFile="mesh.dat";
std::string ibFile="ib.dat";

public:

Interface(Mesh&);
void readInInterface();
void readMesh(Mesh&);
IBDistStat IBLineDist(Node nd, int idx);
wp nearestIBLineDist(Node nd);
void computeLevelSet();
void initializeVolumeFractionField();

};
