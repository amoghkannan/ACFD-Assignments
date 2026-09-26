#pragma once
#include"mesh.h"
#include"vars.h"
#include"mesh.h"
#include"scheme.h"
#include"LinearSolvers/Jacobi.h"
#include"LinearSolvers/GaussSeidel.h"
#include"LinearSolvers/SGS.h"
#include"LinearSolvers/ILU.h"
#include"LinearSolvers/Cholesky.h"
#include"LinearSolvers/GMRES.h"

class Solver{

protected:

public:

Mesh* mesh=nullptr;

Scheme scheme;

//For multigrid
Solver *coarser=nullptr;
Solver *finer=nullptr;

dictionary<std::string,Grid<wp>> vars;

int nVars=0;

dictionary<std::string,Grid<wp>> varsDot;

LinearSolver *ls=nullptr;

Solver();

void setMesh(Mesh& meshIn);
void setMesh(Box box_);
void setVar(std::string name, std::array<idtype,2>indexType_);
void setVar(int ind, Grid<wp>& VarIn);
void setScheme(Scheme& schemeIn);
void setScheme(schemeKey key, schemeVal val);
void setBC(std::string boundary, BCType type, wp val);

Mesh* getMesh();
Scheme& getScheme();

virtual ~Solver(){
        logger.log("Debug: Destroying solver",1);
};

virtual void initialCondition()=0;

virtual void applyBC()=0;

virtual void QDot()=0;

virtual void lSolve(){logger.log("Warning, using undefined lsolve function",1);}; //In case a linear solver is used

virtual void computeTimeStep(Grid<wp>& dt)=0;

virtual void updateVars(Grid<wp>& dt, wp storeFactor)=0;

virtual void updateVars(Grid<wp>& dt, wp storeFactor, std::vector<Grid<wp>>& uStore)=0;

virtual void updateVars(Grid<wp>& dt, wp storeFactor, std::vector<Grid<wp>>& uStore, std::vector<Grid<wp>>& rStore)=0;

virtual wp getResNorm()=0;

virtual bool isConverged()=0;

//For multigrid
virtual Solver* restriction(){ logger.log("Warning, using undefined restriction function",1);return nullptr;};
virtual void prolongation(){ logger.log("Warning, using undefined prolongation function",1);};

};
