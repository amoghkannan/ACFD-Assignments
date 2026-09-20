#pragma once
#include"LinearSolver.h"

class GaussSeidel: public LinearSolver{

protected:

wp rel;

public:

GaussSeidel(int maxIters_, wp tol_, wp rel_, Grid<wp>& varsIn): LinearSolver(maxIters_,tol_,rel_){
        logger.log("Setting up GaussSeidel solver",1);
};
void doIteration(Grid<wp>&) override;
void solve(Grid<wp>&) override;
wp residualCalc(Grid<wp>&) override;
bool isConverged(Grid<wp>&) override;

void invertA(Grid<wp>& var);
};
