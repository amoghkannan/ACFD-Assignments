#pragma once
#include"mesh.h"
#include"vars.h"
#include"mesh.h"
#include"scheme.h"

class Solver{

protected:

Mesh* mesh=nullptr;

std::vector<Grid<wp>> vars;

Scheme scheme;

public:

int nVars=0;

std::vector<Grid<wp>> varsDot;

Solver();

void setMesh(Mesh& meshIn);
void setMesh(int imx, int jmx, int bufW, int bufE, int bufS, int bufN);
std::array<int,6>size();
void setVar();
void setVar(int ind, Grid<wp>& VarIn);
void setScheme(Scheme& schemeIn);
void setScheme(schemeKey key, schemeVal val);
void setBC(int ind, BCType type, wp val);

Mesh* getMesh();
Grid<wp>& getVar(int ind);
Scheme& getScheme();

virtual ~Solver(){
        logger.log("Debug: Destroying solver",1);
};

virtual void initialCondition()=0;

virtual void applyBC()=0;

virtual void QDot()=0;

virtual void computeTimeStep(Grid<wp>& dt)=0;

virtual void updateVars(Grid<wp>& dt, wp storeFactor)=0;

virtual void updateVars(Grid<wp>& dt, wp storeFactor, std::vector<Grid<wp>>& uStore)=0;

virtual void updateVars(Grid<wp>& dt, wp storeFactor, std::vector<Grid<wp>>& uStore, std::vector<Grid<wp>>& rStore)=0;

virtual wp getResNorm()=0;

virtual bool isConverged()=0;

};
