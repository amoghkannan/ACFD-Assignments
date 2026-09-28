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
std::vector<Grid<wp>>res_store;
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
virtual void getResidual()=0;

virtual bool isConverged()=0;

//For multigrid
virtual void setupMultigrid(int level, int maxLevel){logger.log("Invalid multigrid setup",1);};
Grid<wp> restrictVar(Mesh& fineMesh, Grid<wp>& fineVar);
wp volumeAverage(Mesh&fineMesh,Grid<wp>&fineVar,int i,int j);
void prolongateVar(Mesh& fineMesh, Mesh& coarseMesh, Grid<wp>&fineVar, Grid<wp>& coarseVar);
void restriction(std::vector<Grid<wp>>&R_store);
void prolongation(std::vector<Grid<wp>>&R_store);

static wp bilinearInterp2(std::array<Node,2>nds,std::array<wp,2>vals,Node nd);
static wp bilinearInterp4(std::array<Node,4>nds,std::array<wp,4>vals,Node nd);
static Node cc(int i, int j, Mesh& mesh);

};
