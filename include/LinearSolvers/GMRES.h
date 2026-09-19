#include"LinearSolver.h"

class GMRES: public LinearSolver{

protected:

int subspaceSize;
Grid<wp>residualVec;
std::vector<Grid<wp>>krylovVectors;
Grid<wp>H;
std::vector<wp>minCoeffs;
std::vector<wp>minRHS;
wp beta;
Grid<wp>temp;

public:

GMRES(int maxIters_, wp tol_, int ss, Grid<wp>& varsIn): LinearSolver(maxIters_,tol_), subspaceSize(ss), residualVec(varsIn){
        logger.log("Setting up GMRES solver",1);
        for(int i=0;i<=subspaceSize;i++){
                krylovVectors.emplace_back(Grid<wp>(varsIn));
        };

        H=Grid<wp>(subspaceSize+1,subspaceSize,0,0,0,0);
        H.initVal(0.0);
        minCoeffs.resize(subspaceSize);
        minRHS.resize(subspaceSize+1);
};

void doIteration(Grid<wp>&) override;
void solve(Grid<wp>&) override;
wp residualCalc(Grid<wp>&) override;
bool isConverged(Grid<wp>&) override;

void computeResidual(Grid<wp>& var);
void ATimesVec(Grid<wp>& vec, Grid<wp>& ans);
void GrammSchmidt(int idx); 
void givensRotation(); 
void solveForCoeffs();
void updateSoln(Grid<wp>& var);

};
