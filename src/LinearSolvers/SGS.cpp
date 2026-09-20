#include"LinearSolvers/SGS.h"

void SGS::doIteration(Grid<wp>& var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        wp RHS;
        Point p;
        wp data,newElem;

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

                        };
                        
                        var(i,j)=newElem;
                };
        };

        invertAForward(var);

        for(int j=jmx;j>=1;j--){
                for(int i=imx;i>=1;i--){
                        RHS=getRHS(i,j,var);
                        dependencies=getA(i,j,var);
                        newElem=RHS;

                        for(boundMatEntry item:dependencies){
                               p=item.first; 
                               data=item.second;
                              
                               if(p.second<j || (p.second==j && p.first<i)){
                                        newElem=newElem-data*var(p.first,p.second); 
                               };

                        };
                        
                        var(i,j)=newElem;
                };
        };

        invertABackward(var);

};

void SGS::solve(Grid<wp>& var){
        currIter=0;

        logger.log("SGS iteration start: "+std::to_string(residualCalc(var))+"\n",1);

        while(!isConverged(var)){
                logger.log("SGS iteration no: "+std::to_string(currIter)+"\n",1);
                logger.log("SGS residual: "+std::to_string(res)+"\n",1);
                doIteration(var);
                currIter=currIter+1;
        };

};

wp SGS::residualCalc(Grid<wp>& var){
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

bool SGS::isConverged(Grid<wp>& var){
        if(currIter>maxIters || residualCalc(var)<tol) return true;
        return false;
};

void SGS::invertAForward(Grid<wp>& var){
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
                                        continue;
                               };

                        };
                        
                        var(i,j)=var(i,j)/diag;
                };
        };


};

void SGS::invertABackward(Grid<wp>& var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        Point p;
        wp data,diag;

        for(int j=jmx;j>=1;j--){
                for(int i=imx;i>=1;i--){
                        dependencies=getA(i,j,var);

                        for(boundMatEntry item:dependencies){
                               p=item.first; 
                               data=item.second;
                               

                               if(p.second>j || (p.second==j && p.first>i)){
                                        var(i,j)=var(i,j)-data*var(p.first,p.second); 
                               };

                               if(p.first==i && p.second==j){
                                        diag=data;
                                        continue;
                               };

                        };
                        
                        var(i,j)=var(i,j)/diag;
                };
        };


};
