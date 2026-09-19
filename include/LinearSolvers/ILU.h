#include"LinearSolver.h"

class ILU: public LinearSolver{

protected:

wp rel;

public:

ILU(int maxIters_, wp tol_, wp rel_, Grid<wp>& varsIn): LinearSolver(maxIters_,tol_), rel(rel_){
        logger.log("Setting up ILU solver",1);
};
void doIteration(Grid<wp>&) override;
void solve(Grid<wp>&) override;
wp residualCalc(Grid<wp>&) override;
bool isConverged(Grid<wp>&) override;

};
