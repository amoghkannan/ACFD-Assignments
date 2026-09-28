#include"PoissonSolver.h"

std::vector<boundMatEntry> PoissonSolver::getAPoisson(int i,int j,Grid<wp>&var){

        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> ans;
        boundMatEntry temp;

        if(i>1 && i<imx && j>1 && j<jmx){
                temp.first={i-1,j};
                temp.second=-1.0;
                ans.push_back(temp);
        
                temp.first={i+1,j};
                temp.second=-1.0;
                ans.push_back(temp);
        
                temp.first={i,j-1};
                temp.second=-1.0;
                ans.push_back(temp);
        
                temp.first={i,j+1};
                temp.second=-1.0;
                ans.push_back(temp); 
        
                temp.first={i,j};
                temp.second=4.0;
                ans.push_back(temp);

        }
        else{
                temp.first={i,j};
                temp.second=1.0;
                ans.push_back(temp);
        };

        return ans;
};

wp PoissonSolver::getRHSPoisson(int i, int j,Grid<wp>&var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        if(i>1 && i<imx && j>1 && j<jmx){
                return -2.0;
        }
        else{
                return 1.0;
        };
};

PoissonSolver::PoissonSolver(int nx, int ny, int level){
       mesh= new Mesh(Box(nx,ny,1,1,1,1),{NODE,NODE});
       for(int j=1;j<=ny;j++){
                for(int i=1;i<=nx;i++){
                        (*mesh)(i,j)=Node(i*1.0/nx,j*1.0/ny);
                };
       };
       mesh->setGhostNodes();
       
       setVar("phi",{NODE,NODE});
       vars["phi"].initVal(0.0);
       ls=new GaussSeidel(10,1E-4,1.8,vars["phi"]);
       ls->setMatFunc(getAPoisson);
       if(level==1){
                for(int j=1;j<=mesh->intVect[1];j++){
                        for(int i=1;i<=mesh->intVect[0];i++){
                                varsDot["phi"](i,j)=getRHSPoisson(i,j,vars["phi"]);
                        };
                };
                ls->setRHSFunc(varsDot["phi"]);
       };
       this->setScheme(timeStepping,LINEARSOLVER);
       logger.log("Debug: Poisson solver initialization",1);
};

void PoissonSolver::initialCondition(){

        vars["phi"].initVal(0.0);
        logger.log("Debug: Poisson solver initial condition",1);
};

void PoissonSolver::applyBC(){

};

void PoissonSolver::QDot(){

};

void PoissonSolver::computeTimeStep(Grid<wp>& dt){

};

void PoissonSolver::updateVars(Grid<wp>& dt, wp storeFactor){

};

void PoissonSolver::updateVars(Grid<wp>& dt, wp storeFactor, std::vector<Grid<wp>>& uStore){

};

void PoissonSolver::updateVars(Grid<wp>& dt, wp storeFactor, std::vector<Grid<wp>>& uStore, std::vector<Grid<wp>>& rStore){

};

wp PoissonSolver::getResNorm(){
        wp ans=0.0;
        wp RHSCurr;
        Point p;
        wp data,diag,newElem;

        std::vector<boundMatEntry> dependencies;
        
        for(int j=1;j<=mesh->intVect[1];j++){
                for(int i=1;i<=mesh->intVect[0];i++){
                        RHSCurr=varsDot["phi"](i,j);
                        dependencies=getAPoisson(i,j,vars["phi"]);
                        newElem=RHSCurr;

                        for(boundMatEntry item:dependencies){
                               p=item.first; 
                               data=item.second;
                               
                               newElem=newElem-data*vars["phi"](p.first,p.second); 
                        };

                        ans=std::max(ans,fabs(newElem));
                };
        };

        return ans;
};

void PoissonSolver::getResidual(){
        wp RHSCurr;
        Point p;
        wp data,diag,newElem;

        std::vector<boundMatEntry> dependencies;
        
        for(int j=1;j<=mesh->intVect[1];j++){
                for(int i=1;i<=mesh->intVect[0];i++){
                        RHSCurr=varsDot["phi"](i,j);
                        dependencies=getAPoisson(i,j,vars["phi"]);
                        newElem=RHSCurr;

                        for(boundMatEntry item:dependencies){
                               p=item.first; 
                               data=item.second;
                               
                               newElem=newElem-data*vars["phi"](p.first,p.second); 
                        };

                        varsDot["phi"](i,j)=newElem;
                };
        };

};

bool PoissonSolver::isConverged(){
        if(getResNorm()<=1E-4) return true;
        return false;
};

void PoissonSolver::lSolve(){
        
        ls->solve(vars["phi"]);
};

void PoissonSolver::setupMultigrid(int level,int maxLevel){
        if(level==maxLevel) return;
        
        coarser=new PoissonSolver((mesh->box.imx+1)/2,(mesh->box.jmx+1)/2,2);
        coarser->finer=this;
        coarser->mesh=new Mesh((this->mesh)->coarsen());
        coarser->mesh->setGhostNodes();
        coarser->scheme=this->scheme;
        coarser->ls=new GaussSeidel(10,1E-4,1.8,coarser->vars["phi"]);
        (coarser->ls)->setMatFunc(getAPoisson);
        logger.log("Debug: multigrid level set up",1);
        coarser->setupMultigrid(level+1,maxLevel);
};
