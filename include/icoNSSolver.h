#include "solver.h"

class icoNSSolver:public Solver{

private:

std::string meshFile="mesh.dat";
std::string settingsFile="system.dat";
std::string initFile="initial.dat";

public:

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

};
