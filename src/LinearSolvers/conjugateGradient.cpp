#include"LinearSolvers/conjugateGradient.h"

void conjugateGradient::computeResidual(Grid<wp>& var){
        std::array<int,6>sizeArr=residualVec.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        wp RHSCurr;
        Point p;
        wp data,diag,newElem;

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        RHSCurr=getRHS(i,j);
                        dependencies=getA(i,j);
                        newElem=RHSCurr;

                        for(boundMatEntry item:dependencies){
                               p=item.first; 
                               data=item.second;
                               
                               newElem=newElem-data*var(p.first,p.second); 
                        };

                        residualVec(i,j)=newElem; 
                };
        };
       
       if(pp==LEFT_PRECONDITIONER && preconditioner!=nullptr){
                preconditioner->invertA(residualVec);
       };

       if(pp==SPLIT_PRECONDITIONER && preconditioner!=nullptr){
                preconditioner->invertAForward(residualVec);
       };

};

void conjugateGradient::ATimesVec(Grid<wp>& vec, Grid<wp>& ans){
        std::array<int,6>sizeArr=vec.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        Point p;
        wp data,newElem;

       if(pp==RIGHT_PRECONDITIONER && preconditioner!=nullptr){
                preconditioner->invertA(vec);
       };
       if(pp==SPLIT_PRECONDITIONER && preconditioner!=nullptr){
                preconditioner->invertABackward(vec);
       };

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        dependencies=getA(i,j);
                        ans(i,j)=0.0;

                        for(boundMatEntry item:dependencies){
                               p=item.first; 
                               data=item.second;
                               ans(i,j)=ans(i,j)+data*vec(p.first,p.second);
                        };

                };
        };

       if(pp==LEFT_PRECONDITIONER && preconditioner!=nullptr){
                preconditioner->invertA(ans);
       };
       if(pp==SPLIT_PRECONDITIONER && preconditioner!=nullptr){
                preconditioner->invertAForward(ans);
       };


};

void conjugateGradient::doIteration(Grid<wp>& var){

        wp beta;

        rhoNew=pow(residualVec.norm2(),2.0);

        if(currIter==0){
                searchDir=residualVec;
        }
        else{
                beta=rhoNew/rhoOld;
                temp=residualVec*beta;
                searchDir=residualVec+temp;
        };

        ATimesVec(searchDir,temp);
        wp alpha;
        alpha=rhoNew/searchDir.dotProduct(temp);

        temp1=searchDir*alpha;
        var+=temp1;
        temp1=temp*alpha;
        residualVec-=temp1;
        rhoOld=rhoNew;
};

void conjugateGradient::solve(Grid<wp>& var){
        currIter=0;

        std::cout<<"conjugateGradient start; resnorm: "<<std::scientific<<residualCalc(var)<<std::endl;
        
        computeResidual(var);

        while(!isConverged(var)){
                doIteration(var);
                currIter=currIter+1;
                std::cout<<"conjugateGradient step; resnorm: "<<std::scientific<<residualCalc(var)<<std::endl;
        };

        std::cout<<"conjugateGradient iteration finished: "<<std::scientific<<residualCalc(var)<<std::endl;
};

wp conjugateGradient::residualCalc(Grid<wp>& var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        wp RHSCurr;
        Point p;
        wp data,diag,newElem;
        wp ans=0.0;

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        RHSCurr=getRHS(i,j);
                        dependencies=getA(i,j);
                        newElem=RHSCurr;

                        for(boundMatEntry item:dependencies){
                               p=item.first; 
                               data=item.second;
                               
                               newElem=newElem-data*var(p.first,p.second); 
                        };

                        ans=ans+pow(newElem,2.0); 
                };
        };
        
        ans=sqrt(ans);
        res=ans;
        return ans;

};

bool conjugateGradient::isConverged(Grid<wp>& var){
        if(currIter>maxIters || res<tol) return true;
        return false;
};

void conjugateGradient::setupPreconditioner(solverType st, preconditionerPos pp_){
       if(st==JACOBI_SOLVER){
                preconditioner = new Jacobi(maxIters,tol,rel,residualVec,iLim,jLim);
                if(useAExternal){
                         preconditioner->setMatFunc(getAExternal);
                }
                else{
                         preconditioner->setMatFunc(*A);
                };

                if(useRHSExternal){
                         preconditioner->setRHSFunc(getRHSExternal);
                }
                else{
                         preconditioner->setRHSFunc(*RHS);
                };
       }
       else if(st==GAUSS_SEIDEL_SOLVER){
                preconditioner = new GaussSeidel(maxIters,tol,rel,residualVec,iLim,jLim);
                if(useAExternal){
                         preconditioner->setMatFunc(getAExternal);
                }
                else{
                         preconditioner->setMatFunc(*A);
                };

                if(useRHSExternal){
                         preconditioner->setRHSFunc(getRHSExternal);
                }
                else{
                         preconditioner->setRHSFunc(*RHS);
                };

       }
       else if(st==SGS_SOLVER){
                preconditioner = new SGS(maxIters,tol,rel,residualVec,iLim,jLim);
                if(useAExternal){
                         preconditioner->setMatFunc(getAExternal);
                }
                else{
                         preconditioner->setMatFunc(*A);
                };

                if(useRHSExternal){
                         preconditioner->setRHSFunc(getRHSExternal);
                }
                else{
                         preconditioner->setRHSFunc(*RHS);
                };

       }
       else if(st==ILU_SOLVER){
                preconditioner = new ILU(maxIters,tol,rel,residualVec,0,iLim,jLim);
                if(useAExternal){
                         preconditioner->setMatFunc(getAExternal);
                }
                else{
                         preconditioner->setMatFunc(*A);
                };

                if(useRHSExternal){
                         preconditioner->setRHSFunc(getRHSExternal);
                }
                else{
                         preconditioner->setRHSFunc(*RHS);
                };

                (static_cast<ILU*>(preconditioner))->LUDecompose(residualVec);
       }
       else if(st==CHOLESKY_SOLVER){
                preconditioner = new Cholesky(maxIters,tol,rel,residualVec,0,iLim,jLim);
                if(useAExternal){
                         preconditioner->setMatFunc(getAExternal);
                }
                else{
                         preconditioner->setMatFunc(*A);
                };

                if(useRHSExternal){
                         preconditioner->setRHSFunc(getRHSExternal);
                }
                else{
                         preconditioner->setRHSFunc(*RHS);
                };

                (static_cast<Cholesky*>(preconditioner))->LUDecompose(residualVec);
       };

       pp=pp_;

};


