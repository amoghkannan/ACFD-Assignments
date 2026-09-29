#include"icoNSSolver.h"

icoNSSolver::icoNSSolver(){

        std::ifstream infile(meshFile);
        int imx, jmx;

        infile>>imx>>jmx;

        mesh = new Mesh(Box(),{NODE,NODE});

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        infile>>(*mesh)(i,j).x>>(*mesh)(i,j).y;
                };
        };

        mesh->setGhostNodes();
        mesh->calcAreas();
        mesh->calcVolumes();

        infile.close();

        std::ifstream infile1(settingsFile);

        infile1.close();

        (*this).setVar("u",{CELL,CELL});
        (*this).setVar("v",{CELL,CELL});
        (*this).setVar("p",{CELL,CELL});

};

void icoNSSolver::initialCondition(){

        std::ifstream infile(initFile);
        std::string temp;
        infile>>temp;
        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        infile>>vars["u"](i,j);
                };
        };

        infile>>temp;
        for(int j=1;j<=vars["v"].intVect[1];j++){
                for(int i=1;i<=vars["v"].intVect[0];i++){
                        infile>>vars["v"](i,j);
                };
        };

        infile>>temp;
        for(int j=1;j<=vars["p"].intVect[1];j++){
                for(int i=1;i<=vars["p"].intVect[0];i++){
                        infile>>vars["p"](i,j);
                };
        };

        infile.close();
};

void icoNSSolver::applyBC(){

};

void icoNSSolver::QDot(){

};

void icoNSSolver::lSolve(){


};

void icoNSSolver::computeTimeStep(Grid<wp>&dt){

};

void icoNSSolver::updateVars(Grid<wp>&dt, wp storeFactor){

};

void icoNSSolver::updateVars(Grid<wp>&dt, wp storeFactor, std::vector<Grid<wp>>& uStore){

};

void icoNSSolver::updateVars(Grid<wp>&dt, wp storeFactor, std::vector<Grid<wp>>& uStore, std::vector<Grid<wp>>& rStore){

};

wp icoNSSolver::getResNorm(){

};

void icoNSSolver::getResidual(){

};

bool icoNSSolver::isConverged(){

};


