#include"LinearSolvers/Jacobi.h"

void Jacobi::doIteration(Grid<wp>& var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        wp RHS;
        Point p;
        wp data,diag,newElem;

        for(int j=2;j<=jmx-1;j++){
                for(int i=2;i<=imx-1;i++){
                        RHS=getRHS(i,j,var);
                        dependencies=getA(i,j,var);
                        newElem=RHS;

                        for(boundMatEntry item:dependencies){
                               p=item.first; 
                               data=item.second;
                               
                               if(p.first==i && p.second==j){
                                        diag=data;
                                        continue;
                               };

                               newElem=newElem-data*varOld(p.first,p.second); 
                        };

                        var(i,j)=newElem/diag; 
                };
        };
        
        varOld=var;
        applyBC(var);
};

void Jacobi::solve(Grid<wp>& var){
        currIter=0;

        logger.log("Jacobi iteration start: "+std::to_string(residualCalc(var))+"\n",1);
        applyBC(var);

        while(!isConverged(var)){
                logger.log("Jacobi iteration no: "+std::to_string(currIter)+"\n",1);
                logger.log("Jacobi residual: "+std::to_string(res)+"\n",1);
                doIteration(var);
                currIter=currIter+1;
        };

};

wp Jacobi::residualCalc(Grid<wp>& var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        wp RHS;
        Point p;
        wp data,diag,newElem;
        wp ans=0.0;

        for(int j=2;j<=jmx-1;j++){
                for(int i=2;i<=imx-1;i++){
                        RHS=getRHS(i,j,var);
                        dependencies=getA(i,j,var);
                        newElem=RHS;

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

bool Jacobi::isConverged(Grid<wp>& var){
        if(currIter>maxIters || residualCalc(var)<tol) return true;
        return false;
};
