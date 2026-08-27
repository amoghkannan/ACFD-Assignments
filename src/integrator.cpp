#include"integrator.h"

void Integrator::setCFL(wp CFLIn){
        CFL=CFLIn;
};

wp Integrator::getCFL(){
        return CFL;
};

void Integrator::addSolver(Solver& solverIn){

        Solver* newSolver=&solverIn;
        solvers.push_back(newSolver);
 
        std::array<int,6> dims=newSolver->size();

        if(newSolver->getScheme().getScheme(timeStepping)==RK4){
                int n=uStore.size();
                
                if(n==0){
                        for(int i=0;i<newSolver->nVars;i++){
                                uStore.emplace_back(Grid(dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]));
                                rStore.emplace_back(Grid(dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]));
                        };
                }
                else{
                        std::array<int,6> dimsOG=uStore[0].size();
               
                        if(dimsOG[0]<dims[0] || dimsOG[1]<dims[1] || dimsOG[2]<dims[2] ||
                           dimsOG[3]<dims[3] || dimsOG[4]<dims[4] || dimsOG[5]<dims[5]){

                                uStore.clear();
                                rStore.clear();

                                n=n>newSolver->nVars?n:newSolver->nVars;

                                for(int i=0;i<n;i++){
                                        uStore.emplace_back(Grid(dims[0],dims[1],dims[2],
                                                                 dims[3],dims[4],dims[5]));
                                        rStore.emplace_back(Grid(dims[0],dims[1],dims[2],
                                                                 dims[3],dims[4],dims[5]));

                                };
                       
                        }
                        else{
                                if(n<newSolver->nVars){
                                        for(int i=0;i<newSolver->nVars-n;i++){
                                                uStore.emplace_back(Grid(dimsOG[0],dimsOG[1],dimsOG[2],
                                                                        dimsOG[3],dimsOG[4],dimsOG[5]));
                                                rStore.emplace_back(Grid(dimsOG[0],dimsOG[1],dimsOG[2],
                                                                        dimsOG[3],dimsOG[4],dimsOG[5]));

                                        };
                                };
                        };

                };
        
                
        };

        deltaT.emplace_back(Grid(dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]));

        converged.push_back(false);

        nSolvers=nSolvers+1;
};

void Integrator::takeTimeStep(){

       schemeVal timeSteppingScheme;

       Solver* solver;

       for(int i=0;i<nSolvers;i++){

               if(!converged[i]){

               solver=solvers[i];

               timeSteppingScheme=solver->getScheme().getScheme(timeStepping);

               if(timeSteppingScheme==EULER){
                                euler(solver);
               }
               else if(timeSteppingScheme==RK4){
                                rk4(solver);
               }
               else{
                                std::cout<<"Invalid time stepping scheme, exiting"<<std::endl;
                                std::exit(-1);
               };

               if(solver->isConverged()) converged[i]=true;
               };
       };

        nSteps=nSteps+1;

};

void Integrator::integrate(){

       for(int i=0;i<nSolvers;i++){
                solvers[i].initialCondition();
       };

        while(nSteps<maxSteps){
               takeTimeStep(); 
        };
};

void Integrator::rk4(int ind){

        Solver *s=solvers[ind];

        s->computeTimeStep(deltaT[ind]);

        s->QDot();
        s->updateVars(deltaT[ind]);

};

void Integrator::euler(int ind){

        Solver *s=solvers[ind];

        s->computeTimeStep(deltaT[ind]);
        s->QDot();
        s->updateVars(deltaT[ind]);

};
