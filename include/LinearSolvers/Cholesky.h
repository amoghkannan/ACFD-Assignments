#pragma once
#include"LinearSolver.h"
#include<algorithm>

class Cholesky: public LinearSolver{

protected:

wp rel;
int level=0;
Grid<std::vector<boundMatEntry>>LU;
Grid<std::vector<std::pair<int,int>>>S;

public:

Cholesky(int maxIters_, wp tol_, wp rel_, Grid<wp>& varIn, int level_): LinearSolver(maxIters_,tol_,rel_), level(level_){
        logger.log("Setting up Cholesky solver",1);
};
void setLevel(int l){level=l;};
void doIteration(Grid<wp>&) override;
void solve(Grid<wp>&) override;
wp residualCalc(Grid<wp>&) override;
bool isConverged(Grid<wp>&) override;

wp findElement(int iCell,int jCell, int iEntry, int jEntry);
void insertElement(int iCell,int jCell, int iEntry, int jEntry, wp val);
void changeElement(int iCell,int jCell, int iEntry, int jEntry, wp val);
void LUDecompose(Grid<wp>& var);
void invertAForward(Grid<wp>& var) override;
void invertABackward(Grid<wp>& var) override;
};
