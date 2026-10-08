#include"LinearSolver.h"
#include"Jacobi.h"
#include"GaussSeidel.h"
#include"SGS.h"
#include"ILU.h"
#include"Cholesky.h"

class biCGStab: public LinearSolver{

protected:

Grid<wp>temp;
Grid<wp>searchDir;
Grid<wp>residualVec;
Grid<wp>residualVecShadow;
Grid<wp>v;
Grid<wp>s;
Grid<wp>t;

public:

bool terminate=false;
wp rhoOld,rhoNew;
wp alpha;

LinearSolver *preconditioner=nullptr;
preconditionerPos pp;

biCGStab(int maxIters_, wp tol_, wp rel_, Grid<wp>& varsIn, int iLim_, int jLim_): 
        LinearSolver(maxIters_,tol_,rel_, iLim_, jLim_), residualVec(varsIn), temp(varsIn),searchDir(varsIn),
        v(varsIn),s(varsIn),t(varsIn){

        logger.log("Setting up biCGStab solver",1);
};

void doIteration(Grid<wp>&) override;
void solve(Grid<wp>&) override;
wp residualCalc(Grid<wp>&) override;
bool isConverged(Grid<wp>&) override;

void computeResidual(Grid<wp>& var);
void ATimesVec(Grid<wp>& vec, Grid<wp>& ans);

void setupPreconditioner(solverType, preconditionerPos);
};
