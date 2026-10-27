#pragma once
#include"utils.h"
#include"vars.h"
#include"mesh.h"
#include"scheme.h"
#include"solver.h"
#include"icoNSSolver.h"


struct interfaceLine{
        Vec2 l;
        Vec2 n;
        interfaceLine(){};
        interfaceLine(Vec2 l_, Vec2 n_): l(l_),n(n_){};
};

struct interfaceCell{
        int i,j;

        interfaceCell(int i_, int j_): i(i_),j(j_){};
};

struct interfaceFace{
        int bandcellI,bandcellJ;
        int fieldcellI,fieldcellJ;
        int nearestIBLineInd;
        wp nearestIBLineDist;
        Vec2 n;
        char faceType;
        
        interfaceFace(int bI,int bJ,int fI,int fJ,wp nIBD,int nIB,Vec2 n_,char fT):
        bandcellI(bI),bandcellJ(bJ),fieldcellI(fI),fieldcellJ(fJ),nearestIBLineDist(nIBD),nearestIBLineInd(nIB),n(n_),faceType(fT){};

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
std::vector<interfaceCell>bandcells;
std::vector<interfaceCell>interiorcells;
std::vector<interfaceFace>bandfaces;

Grid<Vec2>normals; //Normals to interface
Grid<wp>alphas; //Positions of interface
Grid<wp>F; //Fluxes on I faces
Grid<wp>G; //FLuxes on J faces
Grid<int>cellClass;

Solver *s;
Grid<wp> *IFaceContraVel;
Grid<wp> *JFaceContraVel;

std::string meshFile="mesh.dat";
std::string ibFile="ib.dat";

public:

Interface(Solver&);
void readInInterface();
void readMesh(Mesh&);
IBDistStat IBLineDist(Node nd, int idx);
std::tuple<wp,int,Vec2> nearestIBLineDist(Node nd); //Nearest ibline (signed) dist, index, normal
void computeLevelSet();
void classifyCells();
void IDBandCells();
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
wp enclosedArea(Vec2 m, wp alpha, std::array<Node,4>cellCoords);

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
