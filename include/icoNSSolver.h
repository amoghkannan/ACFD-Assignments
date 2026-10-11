#include "solver.h"

class icoNSSolver:public Solver{

private:

std::string meshFile="mesh.dat";
std::string settingsFile="system.dat";
std::string initFile="init.dat";

dictionary<std::string,wp>params;
wp alpha=0.5;

Grid<wp>tempVar; //Needed for PISO

wp resXMom,resYMom,resMass;

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
void PISODriver();
void computeDiffusiveFlux(std::string varName);
void computeConvectiveFlux();
void computeStoreGrad(Grid<wp>&);
void addPressureSource();
void RhieChow();
void assemblePressureCorrectionEqn();
void correctPressure();
void correctMassFlux();
void correctVelocity();

};
