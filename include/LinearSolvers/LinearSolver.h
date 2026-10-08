#pragma once
#include"../utils.h"
#include"../vars.h"

enum solverType{
        JACOBI_SOLVER,
        GAUSS_SEIDEL_SOLVER,
        SGS_SOLVER,
        ILU_SOLVER,
        CHOLESKY_SOLVER,
        GMRES_SOLVER,
        STEEPEST_DESCENT_SOLVER,
        CONJUGATE_GRADIENT_SOLVER,
        BICONJUGATE_GRADIENT_SOLVER,
        INVALID_SOLVER
};

solverType stringToSolverType(std::string str);

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
int iLim;
int jLim;

public:

boundMatRow getADefault(int i, int j){return (*A)(i,j);};
boundMatRow (*getAExternal)(int,int,int,int)=nullptr;
bool useAExternal=false;

wp getRHSDefault(int i, int j){return (*RHS)(i,j);};
wp (*getRHSExternal)(int,int,int,int)=nullptr;
bool useRHSExternal=false;

Grid<wp>*RHS=nullptr;
Grid<boundMatRow>*A=nullptr;

boundMatRow getA(int i,int j){if(useAExternal) return getAExternal(i,j,iLim,jLim); return getADefault(i,j);};
wp getRHS(int i ,int j){if(useRHSExternal) return getRHSExternal(i,j,iLim,jLim); return getRHSDefault(i,j);};

LinearSolver(int maxIters_, wp tol_, wp rel_, int iLim_, int jLim_): currIter(0), maxIters(maxIters_), tol(tol_), rel(rel_),
                                                                     iLim(iLim_),jLim(jLim_){};
void setMatFunc(boundMatRow(*func)(int,int,int,int)){getAExternal=func;useAExternal=true;logger.log("Setting up A",1);};
void setRHSFunc(wp(*func)(int,int,int,int)){getRHSExternal=func;useRHSExternal=true;logger.log("Setting up RHS",1);};
void setMatFunc(Grid<boundMatRow>&var){A=&var;logger.log("Setting up A",1);};
void setRHSFunc(Grid<wp>&R){RHS=&R;logger.log("Setting up RHS",1);};

virtual void doIteration(Grid<wp>&)=0;
virtual void solve(Grid<wp>&)=0;
virtual wp residualCalc(Grid<wp>&)=0;
virtual bool isConverged(Grid<wp>&)=0;
virtual void invertA(Grid<wp>&){};
virtual void invertAForward(Grid<wp>&){};
virtual void invertABackward(Grid<wp>&){};

};
