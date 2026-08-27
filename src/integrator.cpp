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

        while(nSteps<maxSteps){
               takeTimeStep(); 
        };
};

void Integrator::rk4(Solver* s){


};

void Integrator::euler(Solver* s){


};
