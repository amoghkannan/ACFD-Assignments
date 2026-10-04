#pragma once
#include"solver.h"

class PoissonSolver:public Solver{

private:

public:

static std::vector<boundMatEntry> getAPoisson(int i,int j, int iLim, int jLim);
static wp getRHSPoisson(int i, int j,int iLim, int jLim);

PoissonSolver(int, int,int);

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

void setupMultigrid(int, int) override;

void lSolve() override;

~PoissonSolver() override{
        logger.log("Debug: Destroying Poissonsolver",1);
};

};
