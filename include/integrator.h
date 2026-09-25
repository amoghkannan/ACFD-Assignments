#pragma once
#include"utils.h"
#include"solver.h"
#include<stack>

enum integrationController{
        NO_CONTROLLER,
        MULTIGRID_CONTROLLER
};

class Integrator{

private:

wp CFL=1.0;
integrationController controller;
std::stack<int>turnStack; //Whose turn is it to take a time step?
std::vector<Grid<wp>>deltaT;
std::vector<Grid<wp>>uStore; //For multi-step schemes
std::vector<Grid<wp>>rStore; //For multi-step schemes
dictionary<std::string,Solver*>solvers;
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

void addSolver(Solver& solverIn, std::string name);

void takeTimeStep();
void integrate();

//Time integration schemes

void rk4(int ind);
void euler(int ind);

//Integration controllers
void noController();
void multigridController();

~Integrator(){
        logger.log("Debug: Destroying integrator",1);
};

};
