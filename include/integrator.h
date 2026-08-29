#pragma once
#include"utils.h"
#include"solver.h"
#include<queue>

class Integrator{

private:

wp CFL=1.0;
std::vector<Grid>deltaT;
std::vector<Grid>uStore; //For multi-step schemes
std::vector<Grid>rStore; //For multi-step schemes
std::vector<Solver*>solvers;
std::vector<bool>converged;

int nSolvers=0;
int dumpNumber=0;

void dumpSolution(Solver* s, int ID);
void incDumpNumber();

public:

Integrator();

int nSteps=0;
int maxSteps=1;
int dumpSteps=1;

void setCFL(wp CFLIn);
wp getCFL();

void addSolver(Solver& solverIn);

void takeTimeStep();
void integrate();

//Time integration schemes

void rk4(int ind);
void euler(int ind);

~Integrator(){
        logger.log("Debug: Destroying integrator",1);
};

};
