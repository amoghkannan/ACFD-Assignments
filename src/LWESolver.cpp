#include"LWESolver.h"

LWESolver::LWESolver(Mesh& meshIn){
       setMesh(meshIn);
       setVar();
};

void LWESolver::initialCondition(){


};

void LWESolver::QDot(){


};

void LWESolver::computeTimeStep(Grid& dt){

};

void LWESolver::updateVars(Grid& dt){


};

wp LWESolver::getResNorm(){
        return 0.0;
};

bool LWESolver::isConverged(){

        return true;
};


