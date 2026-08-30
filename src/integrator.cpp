#include"integrator.h"

Integrator::Integrator(){
        logger.log("Debug: Integrator initialization",1);
};

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
        logger.log("Debug: Solver added",1);
};

void Integrator::takeTimeStep(){

       schemeVal timeSteppingScheme;

       Solver *solver;

       for(int i=0;i<nSolvers;i++){

               if(!converged[i]){

               solver=solvers[i];

               timeSteppingScheme=solver->getScheme().getScheme(timeStepping);

               if(timeSteppingScheme==EULER){
                                euler(i);
                                logger.log("Time integration for solver: "+std::to_string(i),LOGCODE);
               }
               else if(timeSteppingScheme==RK4){
                                rk4(i);
                                logger.log("Time integration for solver: "+std::to_string(i),LOGCODE);
               }
               else{
                                std::cout<<"Invalid time stepping scheme, exiting"<<std::endl;
                                std::exit(-1);
               };

               if(solver->isConverged()) converged[i]=true;
               };
       };

        if(nSteps%dumpSteps==0){
               for(int i=0;i<nSolvers;i++){
                        dumpSolution(solvers[i],i);
                };
                incDumpNumber();
        };

        nSteps=nSteps+1;

};

void Integrator::integrate(){

       for(int i=0;i<nSolvers;i++){
                solvers[i]->initialCondition();
                solvers[i]->applyBC();
       };

       logger.log("Debug: Initial condition",1);

        while(nSteps<maxSteps){
               logger.log("Time step number: "+std::to_string(nSteps),LOGCODE);
               takeTimeStep(); 
        };
};

void Integrator::rk4(int ind){

        Solver *s=solvers[ind];

        s->computeTimeStep(deltaT[ind]);
        s->QDot();

        for(int i=0;i<s->nVars;i++){
                uStore[i]=s->getVar(i);
        };

        s->updateVars(deltaT[ind],CFL*0.5);
        
        for(int i=0;i<s->nVars;i++){
                (s->varsDot[i])/6.0;
                rStore[i]=s->varsDot[i];
        };

        s->QDot();

        s->updateVars(deltaT[ind],CFL*0.5,uStore);

        for(int i=0;i<s->nVars;i++){
                (s->varsDot)[i]/3.0;
                rStore[i]+(s->varsDot)[i];
        };

        s->QDot();

        s->updateVars(deltaT[ind],CFL,uStore);

        for(int i=0;i<s->nVars;i++){
                (s->varsDot)[i]/3.0;
                rStore[i]+(s->varsDot)[i];
        };

        s->QDot();

        for(int i=0;i<s->nVars;i++){
                (s->varsDot)[i]/6.0;
                rStore[i]+(s->varsDot)[i];
        };

        s->updateVars(deltaT[ind],CFL,uStore,rStore);

};

void Integrator::euler(int ind){

        Solver *s=solvers[ind];

        s->computeTimeStep(deltaT[ind]);
        s->QDot();
        s->updateVars(deltaT[ind],CFL);
};

void Integrator::dumpSolution(Solver *s, int ID){

        std::ofstream outputfile("s_"+std::to_string(ID)+"_dump_"+std::to_string(dumpNumber)+".dat");
        
        outputfile<<"variables = x y z"<<std::endl;
        int nVars=s->nVars;

        std::array<int,6> dims=s->size();

        int imx=dims[0];
        int jmx=dims[1];

        for(int n=0;n<nVars;n++){
                outputfile<<"var"+std::to_string(n)<<std::endl;       
        };

        outputfile<<"zone T=block0000 i="<<std::to_string(imx)<<" j="<<std::to_string(jmx)<<" k="<<std::to_string(1)<<
                " Datapacking=Block"<<std::endl;
        outputfile<<"Varlocation=([1-"<<std::to_string(nVars+3)<<"]=Nodal)"<<std::endl;
        outputfile<<"STRANDID="<<std::to_string(dumpNumber)<<std::endl;

        for(int j=1;j<=jmx;j++){
           for(int i=1;i<=imx;i++){
                outputfile<<(*(s->getMesh()))(i,j).x<<"\t"<<(*(s->getMesh()))(i,j).y<<"\t"<<0.0<<std::endl;
           }
        };

        for(int n=0;n<nVars;n++){

           for(int j=1;j<=jmx;j++){
              for(int i=1;i<=imx;i++){
                   outputfile<<s->getVar(n)(i,j)<<std::endl;
              }
           };
        };

        outputfile.close();

        std::ofstream logfile("Errors.log",std::ios::app);
        
        while(!logger.isEmpty()){
                logfile<<logger.getLog();
        };
};

void Integrator::incDumpNumber(){
        dumpNumber=dumpNumber+1;
};

