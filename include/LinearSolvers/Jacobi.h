#pragma once
#include"LinearSolver.h"

class Jacobi: public LinearSolver{

protected:

Grid<wp>varOld;

public:

Jacobi(int maxIters_, wp tol_, wp rel_, Grid<wp>& varsIn): LinearSolver(maxIters_,tol_,rel_), varOld(varsIn) {
        logger.log("Setting up Jacobi solver",1);
};
void doIteration(Grid<wp>&) override;
void solve(Grid<wp>&) override;
wp residualCalc(Grid<wp>&) override;
bool isConverged(Grid<wp>&) override;

void invertA(Grid<wp>& var) override;

};
