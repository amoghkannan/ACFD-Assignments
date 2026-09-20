#include"LinearSolvers/GaussSeidel.h"

void GaussSeidel::doIteration(Grid<wp>& var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        wp RHS;
        Point p;
        wp data,newElem;
        wp diag;

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        RHS=getRHS(i,j,var);
                        dependencies=getA(i,j,var);
                        newElem=RHS;

                        for(boundMatEntry item:dependencies){
                               p=item.first; 
                               data=item.second;
                              
                               if(p.second>j || (p.second==j && p.first>i)){
                                        newElem=newElem-data*var(p.first,p.second); 
                               };

                               if(p.first==i && p.second==j){
                                        newElem=newElem+data*var(i,j)*(1.0-rel)/rel;
                               };

                        };
                        
                        var(i,j)=newElem;
                };
        };

        invertA(var);
        
};

void GaussSeidel::solve(Grid<wp>& var){
        currIter=0;

        logger.log("GaussSeidel iteration start: "+std::to_string(residualCalc(var))+"\n",1);

        while(!isConverged(var)){
                logger.log("GaussSeidel iteration no: "+std::to_string(currIter)+"\n",1);
                std::cout<<"GaussSeidel residual: "<<std::scientific<<res<<std::endl;
                doIteration(var);
                currIter=currIter+1;
        };

};

wp GaussSeidel::residualCalc(Grid<wp>& var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        wp RHS;
        Point p;
        wp data,diag,newElem;
        wp ans=0.0;

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
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

bool GaussSeidel::isConverged(Grid<wp>& var){
        if(currIter>maxIters || residualCalc(var)<tol) return true;
        return false;
};

void GaussSeidel::invertA(Grid<wp>& var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        Point p;
        wp data,diag;

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        dependencies=getA(i,j,var);

                        for(boundMatEntry item:dependencies){
                               p=item.first; 
                               data=item.second;
                               

                               if(p.second<j || (p.second==j && p.first<i)){
                                        var(i,j)=var(i,j)-data*var(p.first,p.second); 
                               };

                               if(p.first==i && p.second==j){
                                        diag=data;
                               };

                        };
                        
                        var(i,j)=var(i,j)*rel/diag;
                };
        };

};
