#pragma once
#include"utils.h"
#include"solver.h"
#include<stack>

struct multigridDriver{
};

class Integrator{

private:

wp CFL=1.0;
std::stack<int>turnStack; //Whose turn is it to take a time step?
Grid<wp>deltaT;
std::vector<Grid<wp>>uStore; //For multi-step schemes
std::vector<Grid<wp>>rStore; //For multi-step schemes
dictionary<std::string,Solver*>solvers;
std::vector<bool>converged;

int nSolvers=0;
int dumpNumber=0;

void dumpSolution(Solver* s, int ID);
void incDumpNumber();

int maxMultigridLevel=1;
int timeStepsPerLevel=1;
void multigridController(Solver *s, int level);

public:

Integrator();

int nSteps=0;
int maxSteps=1;
int dumpSteps=1;

void setCFL(wp CFLIn);
wp getCFL();

void addSolver(Solver& solverIn, std::string name);

void takeTimeStep(Solver *solver);
void integrate();

//Time integration schemes

void rk4(Solver *s);
void euler(Solver *s);

~Integrator(){
        logger.log("Debug: Destroying integrator",1);
};

};
