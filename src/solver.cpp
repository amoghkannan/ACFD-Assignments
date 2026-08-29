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

void Solver::setVar(){
       std::array<int,6>dims=this->size();
       vars.emplace_back(Grid(dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]));
       varsDot.emplace_back(Grid(dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]));
       nVars=nVars+1;
};

void Solver::setVar(int ind, Grid& VarIn){
       vars[ind]=VarIn;
};

void Solver::setScheme(Scheme& schemeIn){
        scheme=schemeIn;
};

void Solver::setScheme(schemeKey key, schemeVal val){
        scheme.setScheme(key,val);
};

void Solver::setBC(int ind, BCType type, wp val){
        scheme.setBC(ind,type,val);
};

Mesh* Solver::getMesh(){
        return mesh;
};

Grid& Solver::getVar(int ind){
        if(ind>nVars-1){
                std::cout<<"Invalid variable array access in solver, exiting";
                std::exit(-1);
        };

        return vars[ind];
};

Scheme& Solver::getScheme(){
        return scheme;
};
