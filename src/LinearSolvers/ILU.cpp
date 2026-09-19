#include"LinearSolvers/ILU.h"

void ILU::LUDecompose(Grid<wp>& var){


};

std::vector<boundMatEntry> ILU::getA_ILU(int,int){
        
};

void ILU::doIteration(Grid<wp>& var){
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

                               newElem=newElem-data*var(p.first,p.second); 
                        };

                        var(i,j)=newElem/diag; 
                };
        };
        
        applyBC(var);
};

void ILU::solve(Grid<wp>& var){
        currIter=0;

        logger.log("ILU iteration start: "+std::to_string(residualCalc(var))+"\n",1);
        applyBC(var);

        while(!isConverged(var)){
                logger.log("ILU iteration no: "+std::to_string(currIter)+"\n",1);
                logger.log("ILU residual: "+std::to_string(res)+"\n",1);
                doIteration(var);
                currIter=currIter+1;
        };

};

wp ILU::residualCalc(Grid<wp>& var){
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

bool ILU::isConverged(Grid<wp>& var){
        if(currIter>maxIters || residualCalc(var)<tol) return true;
        return false;
};
