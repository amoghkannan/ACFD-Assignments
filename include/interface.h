#pragma once
#include"utils.h"
#include"vars.h"
#include"mesh.h"
#include"scheme.h"
#include"solver.h"
#include"icoNSSolver.h"

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

class Interface: public Solver{

private:
int nIB;

std::vector<interfaceLine>iblines;
Grid<Vec2>normals; //Normals to interface
Grid<wp>alphas; //Positions of interface
Grid<wp>F; //Fluxes on I faces
Grid<wp>G; //FLuxes on J faces
Grid<int>cellClass;

icoNSSolver *s;
Grid<wp> *IFaceContraVel;
Grid<wp> *JFaceContraVel;

std::string meshFile="mesh.dat";
std::string ibFile="ib.dat";

public:

Interface(icoNSSolver&);
void readInInterface();
void readMesh(Mesh&);
IBDistStat IBLineDist(Node nd, int idx);
wp nearestIBLineDist(Node nd);
void computeLevelSet();
Node InfiniteEdgeIntersection(Vec2 l1, Vec2 l2, bool& intersects);
wp shoelaceArea(std::vector<Node>nds);
bool isLeft(Node nd, Vec2 seg);
std::vector<Node> SutherlandHodgman(std::vector<Vec2>&clipPolygon,std::vector<Node>&subjectPolygon);
//Polygon clipping algorithm used to initialize volume fraction field for VoF; see:
//https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
void initializeVolumeFractionField();
Vec2 calcYoungNormal(int i,int j);
Vec2 calcHtFuncNormal(int i,int j,char dir); 
Vec2 calcNormal(int i,int j); //Uses mixed-Young-Centred method 
wp calcAlpha(int i, int j, Vec2 m); //Piecewise linear interface calc for cell (i,j); returns alpha
void PLIC();
void computeFluxes(wp deltaT);
bool isInside(Vec2 m, wp alpha, Node nd);
wp enclosedArea(Vec2 m, wp alpha, std::vector<Node>cellCoords);

void initialCondition() override;

void applyBC() override;

void QDot() override;

void computeTimeStep(Grid<wp>& dt) override;

void updateVars(Grid<wp>& dt, wp storeFactor) override;

void updateVars(Grid<wp>& dt, wp storeFactor, std::vector<Grid<wp>>& uStore) override;

void updateVars(Grid<wp>& dt, wp storeFactor, std::vector<Grid<wp>>& uStore, std::vector<Grid<wp>>& rStore) override;

wp getResNorm() override;
void getResidual() override;

bool isConverged() override;

};
