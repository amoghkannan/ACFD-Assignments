#include"LinearSolvers/Cholesky.h"

wp Cholesky::findElement(int iCell,int jCell, int iEntry, int jEntry){

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

void Cholesky::insertElement(int iCell,int jCell, int iEntry, int jEntry, wp val){

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

void Cholesky::changeElement(int iCell,int jCell, int iEntry, int jEntry, wp val){

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

void Cholesky::LUDecompose(Grid<wp>& var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        LU=Grid<std::vector<boundMatEntry>>(sizeArr[0],sizeArr[1],sizeArr[2],sizeArr[3],sizeArr[4],sizeArr[5]);
        S=Grid<std::vector<std::pair<int,int>>>(sizeArr[0],sizeArr[1],sizeArr[2],sizeArr[3],sizeArr[4],sizeArr[5]);

        boundMatRow currRow,currRowL;

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        currRow=getA(i,j,var);
                        currRowL.clear();
                        for(auto item:currRow){
                                if((item.first).second<j || ((item.first).second==j && (item.first).first<=i)){
                                        currRowL.push_back(item);
                                        S((item.first).first,(item.first).second).push_back({i,j});
                                };
                        };
                        std::sort(currRowL.begin(),currRowL.end(),[](const boundMatEntry& a, const boundMatEntry& b){
                              return (a.first).second<(b.first).second || 
                                     ((a.first).second==(b.first).second && (a.first).first<(b.first).first);
                        });

                        LU(i,j)=currRowL;
                };
        };


        int ind1,ind2,counter;
        Point p1,p2;
        wp data1,data2;
        wp elem,diag;
        boundMatRow dependencies1,dependencies2;
        
        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                       
                        dependencies1=LU(i,j);
                        counter=0;
                        diag=0.0;

                        for(auto item:dependencies1){
                                if((item.first).first==i && (item.first).second==j) break;

                                elem=item.second;
                                diag=diag+elem*elem;
                                ind2=0;
                              
                                dependencies2=LU((item.first).first,(item.first).second);

                                for(ind1=0;ind1<counter;ind1++){
                                       p1=dependencies1[ind1].first; 
                                       data1=dependencies1[ind1].second;
                                       
                                       data2=0.0;
                                       while(true){
                                                p2=dependencies2[ind2].first; 
                                                data2=dependencies2[ind2].second;
                                                if(p2.second>p1.second || p2.second==p1.second && p2.first> p1.first){
                                                        data2=0.0;
                                                        break;
                                                };
                                                if(p1.first==p2.first && p1.second==p2.second) break;
                                                ind2=ind2+1;
                                       };


                                       elem = elem - data1*data2;   
                                };
                              
                                LU(i,j)[counter].second = elem / dependencies2[dependencies2.size()-1].second;
                                counter=counter+1;
                        };
                        
                        LU(i,j)[dependencies1.size()-1].second=sqrt(LU(i,j)[dependencies1.size()-1].second-diag);

                };
        };

        for(auto item:LU(2,2)){
                std::cout<<(item.first).first<<"\t"<<(item.first).second<<"\t"<<item.second<<std::endl;
        };
};

void Cholesky::doIteration(Grid<wp>& var){
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

void Cholesky::solve(Grid<wp>& var){
        currIter=0;

        logger.log("Cholesky iteration start: "+std::to_string(residualCalc(var))+"\n",1);

        LUDecompose(var);
        doIteration(var);

        logger.log("Cholesky iteration end: "+std::to_string(residualCalc(var))+"\n",1);

};

wp Cholesky::residualCalc(Grid<wp>& var){
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

bool Cholesky::isConverged(Grid<wp>& var){
        if(currIter>maxIters || residualCalc(var)<tol) return true;
        return false;
};

void Cholesky::invertAForward(Grid<wp>& var){
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

                                if(p.first==i && p.second==j){
                                        diag=item.second;
                                        continue;
                                };

                               var(i,j)=var(i,j)-data*var(p.first,p.second); 

                        };
                        
                        var(i,j)=var(i,j)/diag;
                };
        };


};

void Cholesky::invertABackward(Grid<wp>& var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        Point p;
        wp data,diag;

        for(int j=jmx;j>=1;j--){
                for(int i=imx;i>=1;i--){

                        for(auto col:S(i,j)){
                                        if(col.first==i && col.second==j){
                                                diag=findElement(col.first,col.second,i,j);
                                                continue;
                                        };

                                        data=findElement(col.first,col.second,i,j);
                                        var(i,j)=var(i,j)-data*var(col.first,col.second); 
                        };

                        var(i,j)=var(i,j)/diag;
                };
        };


};
