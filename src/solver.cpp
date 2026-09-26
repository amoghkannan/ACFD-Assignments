#include"solver.h"

Solver::Solver(){
        logger.log("Debug: Base solver initialization",1);
};

void Solver::setMesh(Mesh& meshIn){
       mesh=&meshIn; 
};

void Solver::setMesh(Box box_){
        mesh=new Mesh(box_, {NODE,NODE});
};

void Solver::setVar(std::string name, std::array<idtype,2>indexType_){
       vars[name]=Grid<wp>(mesh->box,indexType_);
       varsDot[name]=Grid<wp>(mesh->box,indexType_);
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
