#include "solver.h"

class icoNSSolver:public Solver{

private:

std::string meshFile="mesh.dat";
std::string settingsFile="system.dat";
std::string initFile="init.dat";

dictionary<std::string,wp>params;

public:

dictionary<std::string,std::string>globalBCDict;
dictionary<std::string,Grid<boundMatRow>>coeffs;
dictionary<std::string,Grid<wp>>RHS;
Grid<wp>IFaceMDot;
Grid<wp>JFaceMDot;

icoNSSolver();
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

void SIMPLEDriver();
void computeDiffusiveFlux(std::string varName);
void computeConvectiveFlux();
void RhieChow();

};
