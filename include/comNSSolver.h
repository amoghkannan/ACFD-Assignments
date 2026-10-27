#include "solver.h"

class comNSSolver:public Solver{

private:

std::string meshFile="mesh.dat";
std::string settingsFile="system.dat";
std::string initFile="init.dat";

dictionary<std::string,wp>params;

wp resXMom,resYMom,resMass,resEnergy;

public:

dictionary<std::string,Grid<wp>>x_qp_left;
dictionary<std::string,Grid<wp>>x_qp_right;
dictionary<std::string,Grid<wp>>y_qp_left;
dictionary<std::string,Grid<wp>>y_qp_right;

dictionary<std::string,Grid<wp>>F;
dictionary<std::string,Grid<wp>>G;

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
std::vector<wp>conservativeToPrimitive(std::vector<wp>qc, wp gm, wp RGas);
std::vector<wp>primitiveToConservative(std::vector<wp>qp, wp gm, wp RGas);

};
