#pragma once
#include"../utils.h"
#include"../vars.h"

enum solverType{
        JACOBI_SOLVER,
        GAUSS_SEIDEL_SOLVER,
        SGS_SOLVER,
        ILU_SOLVER,
        CHOLESKY_SOLVER,
        GMRES_SOLVER
};

enum preconditionerPos{
        LEFT_PRECONDITIONER,
        RIGHT_PRECONDITIONER,
        SPLIT_PRECONDITIONER
};

class LinearSolver{

protected:

int currIter;
int maxIters;
wp rel;
wp tol;
wp res;

std::vector<boundMatEntry>(*getA)(int,int,Grid<wp>&)=nullptr;
wp (*getRHS)(int,int,Grid<wp>&)=nullptr;


public:

LinearSolver(int maxIters_, wp tol_, wp rel_): currIter(0), maxIters(maxIters_), tol(tol_), rel(rel_){};
void setMatFunc(std::vector<boundMatEntry>(*func)(int,int,Grid<wp>&)){getA=func;logger.log("Setting up A",1);};
void setRHSFunc(wp (*func)(int,int,Grid<wp>&)){getRHS=func;logger.log("Setting up RHS",1);};

virtual void doIteration(Grid<wp>&)=0;
virtual void solve(Grid<wp>&)=0;
virtual wp residualCalc(Grid<wp>&)=0;
virtual bool isConverged(Grid<wp>&)=0;
virtual void invertA(Grid<wp>&){};
virtual void invertAForward(Grid<wp>&){};
virtual void invertABackward(Grid<wp>&){};

};
