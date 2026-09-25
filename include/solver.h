#pragma once
#include"mesh.h"
#include"vars.h"
#include"mesh.h"
#include"scheme.h"

class Solver{

protected:

Mesh* mesh=nullptr;

Scheme scheme;

public:

dictionary<std::string,Grid<wp>> vars;

int nVars=0;

dictionary<std::string,Grid<wp>> varsDot;

Solver();

void setMesh(Mesh& meshIn);
void setMesh(int imx, int jmx, int bufW, int bufE, int bufS, int bufN);
std::array<int,6>size();
void setVar(std::string name);
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

virtual void computeTimeStep(Grid<wp>& dt)=0;

virtual void updateVars(Grid<wp>& dt, wp storeFactor)=0;

virtual void updateVars(Grid<wp>& dt, wp storeFactor, std::vector<Grid<wp>>& uStore)=0;

virtual void updateVars(Grid<wp>& dt, wp storeFactor, std::vector<Grid<wp>>& uStore, std::vector<Grid<wp>>& rStore)=0;

virtual wp getResNorm()=0;

virtual bool isConverged()=0;

};
