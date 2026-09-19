#include"LinearSolver.h"

class Jacobi: public LinearSolver{

protected:

wp rel;

public:

Jacobi(int maxIters_, wp tol_, wp rel_): LinearSolver(maxIters_,tol_), rel(rel_) {};
void doIteration() override;
void solve(Grid<wp>&) override;
wp residualCalc() override;
bool isConverged() override;

};
