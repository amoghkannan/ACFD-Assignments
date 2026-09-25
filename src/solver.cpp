#include"solver.h"

Solver::Solver(){
        logger.log("Debug: Base solver initialization",1);
};

void Solver::setMesh(Mesh& meshIn){
       mesh=&meshIn; 
};

void Solver::setMesh(int imx, int jmx, int bufW, int bufE, int bufS, int bufN){
        mesh=new Mesh(imx,jmx,bufW,bufE,bufS,bufN);
};

std::array<int,6> Solver::size(){
        return mesh->size();
};

void Solver::setVar(std::string name){
       std::array<int,6>dims=this->size();
       vars[name]=Grid<wp>(dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]);
       varsDot[name]=Grid<wp>(dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]);
       nVars=nVars+1;
};

void Solver::setVar(int ind, Grid<wp>& VarIn){
       vars(ind)=VarIn;
};

void Solver::setScheme(Scheme& schemeIn){
        scheme=schemeIn;
};

void Solver::setScheme(schemeKey key, schemeVal val){
        scheme.setScheme(key,val);
};

void Solver::setBC(std::string boundary, BCType type, wp val){
        scheme.setBC(boundary,type,val);
};

Mesh* Solver::getMesh(){
        return mesh;
};

Scheme& Solver::getScheme(){
        return scheme;
};
