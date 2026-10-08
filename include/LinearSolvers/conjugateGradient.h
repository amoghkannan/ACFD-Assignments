#include"LinearSolver.h"
#include"Jacobi.h"
#include"GaussSeidel.h"
#include"SGS.h"
#include"ILU.h"
#include"Cholesky.h"

class conjugateGradient: public LinearSolver{

protected:

Grid<wp>temp;
Grid<wp>temp1;
Grid<wp>searchDir;
Grid<wp>residualVec;

wp rhoOld,rhoNew;

public:

LinearSolver *preconditioner=nullptr;
preconditionerPos pp;

conjugateGradient(int maxIters_, wp tol_, wp rel_, Grid<wp>& varsIn, int iLim_, int jLim_): 
        LinearSolver(maxIters_,tol_,rel_, iLim_, jLim_), residualVec(varsIn), temp(varsIn){

        logger.log("Setting up conjugateGradient solver",1);
};

void doIteration(Grid<wp>&) override;
void solve(Grid<wp>&) override;
wp residualCalc(Grid<wp>&) override;
bool isConverged(Grid<wp>&) override;

void computeResidual(Grid<wp>& var);
void ATimesVec(Grid<wp>& vec, Grid<wp>& ans);

void setupPreconditioner(solverType, preconditionerPos);
};
