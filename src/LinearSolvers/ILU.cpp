#include"LinearSolvers/ILU.h"

wp ILU::findElement(int iCell,int jCell, int iEntry, int jEntry){

        std::vector<boundMatEntry> dependencies;
        Point p;
        wp data,ans;

        dependencies=LU(iCell,jCell);

        ans=0.0;

        for(boundMatEntry item:dependencies){
               p=item.first; 
               data=item.second;
               
               if(p.first==iEntry && p.second==jEntry){
                         ans=data;
                         break;
               };
        };

        return ans;

};

void ILU::insertElement(int iCell,int jCell, int iEntry, int jEntry, wp val){

        std::vector<boundMatEntry> dependencies;
        Point p;
        wp data;

        dependencies=LU(iCell,jCell);

        bool isFound=false;
        int counter=0;

        for(boundMatEntry item:dependencies){
               p=item.first; 
               data=item.second;
               
               if(p.first==iEntry && p.second==jEntry){
                         isFound=true;
                         LU(iCell,jCell)[counter].second=val;
                         break;
               };

               counter=counter+1;
        };

        if(!isFound && level>0){
                boundMatEntry newEntry={{iEntry,jEntry},val};
                LU(iCell,jCell).push_back(newEntry);
        };
};

void ILU::changeElement(int iCell,int jCell, int iEntry, int jEntry, wp val){

        std::vector<boundMatEntry> dependencies;
        Point p;
        wp data;

        dependencies=LU(iCell,jCell);

        bool isFound=false;
        int counter=0;

        for(boundMatEntry item:dependencies){
               p=item.first; 
               data=item.second;
               
               if(p.first==iEntry && p.second==jEntry){
                         isFound=true;
                         LU(iCell,jCell)[counter].second+=val;
                         break;
               };

               counter=counter+1;
        };

        if(!isFound && level>0){
                boundMatEntry newEntry={{iEntry,jEntry},val};
                LU(iCell,jCell).push_back(newEntry);
        };
};

void ILU::LUDecompose(Grid<wp>& var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        LU=Grid<std::vector<boundMatEntry>>(sizeArr[0],sizeArr[1],sizeArr[2],sizeArr[3],sizeArr[4],sizeArr[5]);
        
        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        LU(i,j)=getA(i,j,var);
                };
        };

        std::vector<boundMatEntry> dependencies;
        Point p;
        wp data,diag;
        wp factor;
        wp elem;


        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        
                        diag=findElement(i,j,i,j);

                        dependencies=LU(i,j);

                        for(int i1=i+1;i1<=imx;i1++){
                                
                                factor=findElement(i1,j,i,j)/diag;

                                if(factor!=0.0){
                                        for(boundMatEntry item:dependencies){
                                               p=item.first; 
                                               data=item.second;

                                               if(p.second>j || (p.second==j && p.first>i)){
                                                      changeElement(i1,j,p.first,p.second,-factor*data);  
                                               };
                                               
                                        };
                                };
                                
                                insertElement(i1,j,i,j,factor);

                        };

                        for(int j1=j+1;j1<=jmx;j1++){
                                for(int i1=1;i1<=imx;i1++){

                                        factor=findElement(i1,j1,i,j)/diag;

                                        if(factor!=0.0){
                                                for(boundMatEntry item:dependencies){
                                                       p=item.first; 
                                                       data=item.second;

                                                       if(p.second>j || (p.second==j && p.first>i)){
                                                              changeElement(i1,j1,p.first,p.second,-factor*data);  
                                                       };
                                                       
                                                };
                                        };
                                        
                                        insertElement(i1,j1,i,j,factor);


                                };
                        };

                };
        };
        
};

void ILU::doIteration(Grid<wp>& var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        wp RHS;

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        RHS=getRHS(i,j,var);
                        var(i,j)=RHS;

                };
        };

        invertAForward(var);
        invertABackward(var);
       
};

void ILU::solve(Grid<wp>& var){
        currIter=0;

        logger.log("ILU iteration start: "+std::to_string(residualCalc(var))+"\n",1);

        LUDecompose(var);
        doIteration(var);

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

bool ILU::isConverged(Grid<wp>& var){
        if(currIter>maxIters || residualCalc(var)<tol) return true;
        return false;
};

void ILU::invertAForward(Grid<wp>& var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        Point p;
        wp data,diag;

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        dependencies=LU(i,j);

                        for(boundMatEntry item:dependencies){
                               p=item.first; 
                               data=item.second;
                               

                               if(p.second<j || (p.second==j && p.first<i)){
                                        var(i,j)=var(i,j)-data*var(p.first,p.second); 
                               };


                        };
                        
                        diag=1.0;
                        var(i,j)=var(i,j)/diag;
                };
        };


};

void ILU::invertABackward(Grid<wp>& var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        Point p;
        wp data,diag;

        for(int j=jmx;j>=1;j--){
                for(int i=imx;i>=1;i--){
                        dependencies=LU(i,j);

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
