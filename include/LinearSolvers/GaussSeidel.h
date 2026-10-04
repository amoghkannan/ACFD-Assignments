#pragma once
#include"LinearSolver.h"

class GaussSeidel: public LinearSolver{

protected:

public:

GaussSeidel(int maxIters_, wp tol_, wp rel_, Grid<wp>& varsIn, int iLim_, int jLim_): 
        LinearSolver(maxIters_,tol_,rel_,iLim_,jLim_){
        logger.log("Setting up GaussSeidel solver",1);
};
void doIteration(Grid<wp>&) override;
void solve(Grid<wp>&) override;
wp residualCalc(Grid<wp>&) override;
bool isConverged(Grid<wp>&) override;

void invertA(Grid<wp>& var) override;
};
