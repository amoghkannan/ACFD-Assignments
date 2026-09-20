#pragma once
#include"LinearSolver.h"

class SGS: public LinearSolver{

protected:

wp rel;

public:

SGS(int maxIters_, wp tol_, wp rel_, Grid<wp>& varsIn): LinearSolver(maxIters_,tol_,rel_){
        logger.log("Setting up SGS solver",1);
};
void doIteration(Grid<wp>&) override;
void solve(Grid<wp>&) override;
wp residualCalc(Grid<wp>&) override;
bool isConverged(Grid<wp>&) override;

void invertAForward(Grid<wp>& var) override;
void invertABackward(Grid<wp>& var) override;

};
