#include "solver.h"
#include"interface.h"

class comNSSolver:public Solver{

private:

std::string meshFile="mesh.dat";
std::string settingsFile="system.dat";
std::string initFile="init.dat";

dictionary<std::string,std::pair<int,int>>lo;
dictionary<std::string,std::pair<int,int>>hi;
dictionary<std::string,std::pair<int,int>>inc;
dictionary<std::string,std::pair<int,int>>inc1;
dictionary<std::string,wp>sig;


wp resXMom,resYMom,resMass,resEnergy;

public:

dictionary<std::string,wp>params;

Interface *in=nullptr;

dictionary<std::string,Grid<wp>>x_qp_left;
dictionary<std::string,Grid<wp>>x_qp_right;
dictionary<std::string,Grid<wp>>y_qp_left;
dictionary<std::string,Grid<wp>>y_qp_right;

dictionary<std::string,Grid<wp>>F;
dictionary<std::string,Grid<wp>>G;

Grid<wp>temperature;
Grid<wp>mu;

dictionary<std::string,wp>res;
dictionary<std::string,wp>resScale;

dictionary<std::string,std::string>globalBCDict;
dictionary<std::string,bool>zeroFluxBC;

comNSSolver();
void initialCondition() override;
void applyBC() override;
void QDot() override;
void lSolve() override;
void computeTimeStep(Grid<wp>&dt) override;
void updateVars(Grid<wp>&dt, wp storeFactor) override;
void updateVars(Grid<wp>&dt, wp storeFactor, std::vector<Grid<wp>>& uStore) override;
void updateVars(Grid<wp>&dt, wp storeFactor, std::vector<Grid<wp>>& uStore, std::vector<Grid<wp>>& rStore) override;
wp getResNorm() override;
void getResidual() override;
bool isConverged() override;

void reconstructStates();
void computeFluxes();
void conservativeToPrimitive(std::vector<wp>&qc, wp gm, wp RGas);
void primitiveToConservative(std::vector<wp>&qp, wp gm, wp RGas);
void computeTemperature(Grid<wp>&pressure,Grid<wp>&density,wp rgas);
void computeViscosity();
};
