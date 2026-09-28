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

void Integrator::addSolver(Solver& solverIn, std::string name){
        Solver* newSolver=&solverIn;
        solvers[name]=newSolver;
 
        if(newSolver->getScheme().getScheme(timeStepping)==RK4 || newSolver->getScheme().schemeDict.hasKey(multigrid)){
                int n=uStore.size();
                if(n<newSolver->nVars){
                        uStore.resize(newSolver->nVars);
                        rStore.resize(newSolver->nVars);
                };

                for(int i=0;i<n;i++){
                        uStore[i]=newSolver->vars(i);
                        rStore[i]=newSolver->vars(i);
                };

        };
        if(newSolver->getScheme().getScheme(timeStepping)==RK4 || newSolver->getScheme().getScheme(timeStepping)==EULER){
                deltaT=newSolver->vars(0);
        };

         if(newSolver->getScheme().schemeDict.hasKey(multigrid)){
                newSolver->setupMultigrid(1,maxMultigridLevel);
        };
 
        converged.push_back(false);

        nSolvers=nSolvers+1;
        logger.log("Debug: Solver added",1);
};

void Integrator::takeTimeStep(Solver *solver){

       schemeVal timeSteppingScheme;

       timeSteppingScheme=solver->getScheme().getScheme(timeStepping);

       if(timeSteppingScheme==EULER){
                        euler(solver);
       }
       else if(timeSteppingScheme==RK4){
                        rk4(solver);
       }
       else if(timeSteppingScheme==LINEARSOLVER){
                        solver->lSolve();
       }
       else{
                        std::cout<<"Invalid time stepping scheme, exiting"<<std::endl;
                        std::exit(-1);
       };

};

void Integrator::integrate(){

        auto start=std::chrono::high_resolution_clock::now();

       for(int i=0;i<nSolvers;i++){
                solvers(i)->initialCondition();
                solvers(i)->applyBC();
       };

       logger.log("Debug: Initial condition",1);

        int notConvergedCounter;

        while(nSteps<maxSteps){
               logger.log("Time step number: "+std::to_string(nSteps),LOGCODE);
               notConvergedCounter=0;
               for(int i=0;i<nSolvers;i++){
                        if(converged[i]) continue;
                        notConvergedCounter++;
                        if(solvers(i)->getScheme().hasScheme(multigrid)){
                                multigridController(solvers(i),1);
                        }
                        else{
                                takeTimeStep(solvers(i));
                        };
                        if(solvers(i)->isConverged()) converged[i]=true;
               };

               nSteps=nSteps+1;
               if(nSteps%dumpSteps==0){
                      for(int i=0;i<nSolvers;i++){
                               dumpSolution(solvers(i),i);
                       };
                       incDumpNumber();
               };
               
               if(notConvergedCounter==0) break;
        };

        logger.log("Finished!",1);

        auto end=std::chrono::high_resolution_clock::now();

        std::chrono::duration<double>elapsed=end-start;

        std::cout<<"Solved finished task in "<<elapsed.count()<<" seconds"<<std::endl;
};

void Integrator::rk4(Solver *s){

        s->computeTimeStep(deltaT);
        s->QDot();

        for(int i=0;i<s->nVars;i++){
                uStore[i]=s->vars(i);
        };

        s->updateVars(deltaT,CFL*0.5);
        
        for(int i=0;i<s->nVars;i++){
                (s->varsDot(i))/=6.0;
                rStore[i]=(s->varsDot)(i);
        };

        s->QDot();

        s->updateVars(deltaT,CFL*0.5,uStore);

        for(int i=0;i<s->nVars;i++){
                (s->varsDot)(i)/=3.0;
                rStore[i]+=(s->varsDot)(i);
        };

        s->QDot();

        s->updateVars(deltaT,CFL,uStore);

        for(int i=0;i<s->nVars;i++){
                (s->varsDot)(i)/=3.0;
                rStore[i]+=(s->varsDot)(i);
        };

        s->QDot();

        for(int i=0;i<s->nVars;i++){
                (s->varsDot)(i)/=6.0;
                rStore[i]+=(s->varsDot)(i);
        };

        s->updateVars(deltaT,CFL,uStore,rStore);

};

void Integrator::euler(Solver *s){

        s->computeTimeStep(deltaT);
        s->QDot();
        s->updateVars(deltaT,CFL);
};

void Integrator::multigridController(Solver *s, int level){
        
        switch(s->getScheme().getScheme(multigrid)){
                case(V_CYCLE):
                        logger.log("Multigrid level: "+std::to_string(level),1);
                        takeTimeStep(s);
                        if(level!=maxMultigridLevel){
                                logger.log("Multigrid level: "+std::to_string(level)+" restriction",1);
                                s->restriction(rStore);
                                multigridController(s->coarser,level+1);
                                logger.log("Multigrid level: "+std::to_string(level)+" prolongation",1);
                                s->coarser->prolongation(rStore);
                                logger.log("Multigrid level: "+std::to_string(level)+" final relaxation",1);
                                takeTimeStep(s);
                        };
                        return;
                case(W_CYCLE):
                        takeTimeStep(s);
                        if(level!=maxMultigridLevel){
                                s->restriction(rStore);
                                multigridController(s->coarser,level+1);
                                s->coarser->prolongation(rStore);
                                takeTimeStep(s);
                                s->restriction(rStore);
                                multigridController(s->coarser,level+1);
                                s->coarser->prolongation(rStore);
                                takeTimeStep(s);
                        };
                        return;
                default:
                        logger.log("Error, invalid multigrid scheme",1);
                        break;

        };
};

void Integrator::dumpSolution(Solver *s, int ID){

        std::ofstream outputfile("s_"+std::to_string(ID)+"_dump_"+std::to_string(dumpNumber)+".dat");
        
        outputfile<<"variables = x y z"<<std::endl;
        int nVars=s->nVars;

        std::array<int,6> dims=s->mesh->size();

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
                   outputfile<<s->vars(n)(i,j)<<std::endl;
              }
           };
        };

        outputfile.close();

        std::ofstream logfile("Errors.log",std::ios::app);
        
        while(!logger.isEmpty()){
                logfile<<logger.getLog()<<std::endl;
        };
};

void Integrator::incDumpNumber(){
        dumpNumber=dumpNumber+1;
};

