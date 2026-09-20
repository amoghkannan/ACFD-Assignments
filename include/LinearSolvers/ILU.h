#include"LinearSolver.h"

class ILU: public LinearSolver{

protected:

wp rel;
int level=0;
std::vector<boundMatRow>LU;

public:

ILU(int maxIters_, wp tol_, wp rel_, Grid<wp>& varIn): LinearSolver(maxIters_,tol_,rel_){
        logger.log("Setting up ILU solver",1);
};
void setLevel(int l){level=l;};
void doIteration(Grid<wp>&) override;
void solve(Grid<wp>&) override;
wp residualCalc(Grid<wp>&) override;
bool isConverged(Grid<wp>&) override;

void LUDecompose(Grid<wp>& var);
std::vector<boundMatEntry>getA_ILU(int,int);
};
