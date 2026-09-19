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
wp res;

void (*applyBC)()=nullptr;
std::vector<boundMatEntry>(*getA)(int,int)=nullptr;
wp (*getRHS)(int,int)=nullptr;

public:

LinearSolver(int maxIters_, wp tol_): currIter(0), maxIters(maxIters_), tol(tol_) {};
void setBCFunc(void (*func)()){applyBC=func;logger.log("Setting up BCFunc",1);};
void setMatFunc(std::vector<boundMatEntry>(*func)(int,int)){getA=func;logger.log("Setting up A",1);};
void setRHSFunc(wp (*func)(int,int)){getRHS=func;logger.log("Setting up RHS",1);};

virtual void doIteration(Grid<wp>&)=0;
virtual void solve(Grid<wp>&)=0;
virtual wp residualCalc(Grid<wp>&)=0;
virtual bool isConverged(Grid<wp>&)=0;
};
