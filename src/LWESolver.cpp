#include"LWESolver.h"

LWESolver::LWESolver(Mesh& meshIn){
       setMesh(meshIn);
       setVar();
       logger.log("Debug: LWE solver initialization",1);
};

void LWESolver::initialCondition(){

       wp deltaX=(*mesh)(2,1).x-(*mesh)(1,1).x;

       std::array<int,6>dims=this->size();

        for(int j=1;j<=dims[1];j++){
                for(int i=1;i<=dims[0];i++){
                        vars[0](i,j)=sin(M_PI*(*mesh)(i,j).x/(10.0*deltaX));
                };
        };

        logger.log("Debug: LWE solver initial condition",1);
};

void LWESolver::applyBC(){
       std::array<int,6>dims=this->size();
       
       int imx=dims[0];
       int jmx=dims[1];

       vars[0](0,1)=vars[0](imx,1);
       vars[0](imx+1,1)=vars[0](1,1);
};

void LWESolver::QDot(){

       schemeVal derivativeXScheme=scheme.getScheme(derivativeX);

       switch(derivativeXScheme){
                case FOU:
                        if(c>0.0){
                                scheme.EBD1(vars[0],varsDot[0],*mesh,'x');
                                varsDot[0]*=(-c);
                        }
                        else{
                                scheme.EFD1(vars[0],varsDot[0],*mesh,'x');
                                varsDot[0]*=(-c);
                        };
                        break;
                case C4:
                        scheme.ICD4(vars[0],varsDot[0],*mesh,'x');
                        varsDot[0]*=(-c);
                        break;
                default:
                        break;

       };
};

void LWESolver::computeTimeStep(Grid<wp>& dt){
 
       wp deltaX=(*mesh)(2,1).x-(*mesh)(1,1).x;

       std::array<int,6>dims=dt.size();

        for(int j=1;j<=dims[1];j++){
                for(int i=1;i<=dims[0];i++){
                        dt(i,j) = deltaX / fabs(c);
                };
        };
       
};

void LWESolver::updateVars(Grid<wp>& dt, wp storeFactor){

       std::array<int,6>dims=dt.size();

        for(int j=1;j<=dims[1];j++){
                for(int i=1;i<=dims[0];i++){
                        vars[0](i,j)= vars[0](i,j)+dt(i,j)*storeFactor*varsDot[0](i,j);
                };
        };

        applyBC();
};

void LWESolver::updateVars(Grid<wp>& dt, wp storeFactor, std::vector<Grid<wp>>& uStore){

       wp deltaX=(*mesh)(2,1).x-(*mesh)(1,1).x;

       std::array<int,6>dims=dt.size();

        for(int j=1;j<=dims[1];j++){
                for(int i=1;i<=dims[0];i++){
                        vars[0](i,j)= uStore[0](i,j)+dt(i,j)*storeFactor*varsDot[0](i,j);
                };
        };
        applyBC();
};

void LWESolver::updateVars(Grid<wp>& dt, wp storeFactor, std::vector<Grid<wp>>& uStore, std::vector<Grid<wp>>& rStore){

       wp deltaX=(*mesh)(2,1).x-(*mesh)(1,1).x;

       std::array<int,6>dims=dt.size();

        for(int j=1;j<=dims[1];j++){
                for(int i=1;i<=dims[0];i++){
                        vars[0](i,j)= uStore[0](i,j)+dt(i,j)*storeFactor*rStore[0](i,j);
                };
        };

        applyBC();

};

wp LWESolver::getResNorm(){
        return 0.0;
};

bool LWESolver::isConverged(){

        return false;
};


