#include"../utils.h"
#include"../vars.h"

enum solverType{
        JACOBI
};


class LinearSolver{

protected:

int currIter;
int maxIters;
wp tol;

void (*applyBC)()=nullptr;
std::vector<boundMatEntry>(*getA)(int,int)=nullptr;
wp (*getRHS)(int,int)=nullptr;

public:

LinearSolver(int maxIters_, wp tol_): currIter(0), maxIters(maxIters_), tol(tol_) {};
void setBCFunc(void (*func)()){applyBC=func;};
void setMatFunc(std::vector<boundMatEntry>(*func)(int,int)){getA=func;};
void setRHSFunc(wp (*func)(int,int)){getRHS=func;};

virtual void doIteration()=0;
virtual void solve(Grid<wp>&)=0;
virtual wp residualCalc()=0;
virtual bool isConverged()=0;
};
