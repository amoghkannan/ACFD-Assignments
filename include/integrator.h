#pragma once
#include"utils.h"
#include"solver.h"
#include<queue>

class Integrator{

private:

wp CFL=1.0;
std::vector<Grid>deltaT;
std::vector<Grid>uStore; //For multi-step schemes
std::vector<Grid>dotStore; //For multi-step schemes
std::vector<Solver*>solvers={nullptr};
std::vector<bool>converged;

int nSolvers=0;

std::queue<int>updateOrder;


public:

int nSteps=0;
int maxSteps=1;

void setCFL(wp CFLIn);
wp getCFL();

void addSolver(Solver& solverIn);

void takeTimeStep();
void integrate();

//Time integration schemes

void rk4(int ind);
void euler(int ind);

};
