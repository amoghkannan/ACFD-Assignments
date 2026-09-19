#include"LinearSolvers/GMRES.h"

void GMRES::computeResidual(Grid<wp>& var){
        std::array<int,6>sizeArr=residualVec.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        wp RHS;
        Point p;
        wp data,diag,newElem;

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

                        residualVec(i,j)=newElem; 
                };
        };
        
};

void GMRES::ATimesVec(Grid<wp>& vec, Grid<wp>& ans){
        std::array<int,6>sizeArr=vec.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

        std::vector<boundMatEntry> dependencies;
        Point p;
        wp data,newElem;

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        dependencies=getA(i,j,vec);
                        ans(i,j)=0.0;

                        for(boundMatEntry item:dependencies){
                               p=item.first; 
                               data=item.second;
                               ans(i,j)=ans(i,j)+data*vec(p.first,p.second);
                        };

                };
        };

};

void GMRES::GrammSchmidt(int idx){
       
        for(int i=0;i<idx;i++){
                H(i+1,idx)=krylovVectors[idx].dotProduct(krylovVectors[i]);
                temp=krylovVectors[i];
                temp*=H(i+1,idx);
                krylovVectors[idx]-=temp;
        };

        H(idx+1,idx)=krylovVectors[idx].norm2();

};

void GMRES::givensRotation(){
       
       wp cos_theta;
       wp sin_theta;
       wp temp1;
       wp temp2;

       for(int i=0;i<=subspaceSize;i++){
                minRHS[i]=0.0;
       };

       minRHS[0]=beta;

       for(int i=1;i<=subspaceSize;i++){
          cos_theta=H(i,i)/sqrt(pow(H(i,i),2)+pow(H(i+1,i),2));
          sin_theta=-H(i+1,i)/sqrt(pow(H(i,i),2)+pow(H(i+1,i),2));

          for(int j=i;j<=subspaceSize;j++){
                temp1=H(i,j);
                temp2=H(i+1,j);

                H(i,j)=temp1*cos_theta-temp2*sin_theta;
                H(i+1,j)=temp1*sin_theta+temp2*cos_theta;
          };

          temp1=minRHS[i-1];
          minRHS[i-1]=temp1*cos_theta;
          minRHS[i]=temp1*sin_theta;
       };

};

void GMRES::solveForCoeffs(){

        for(int i=subspaceSize;i>=1;i--){
               minCoeffs[i-1]=minRHS[i-1]; 
               for(int j=i+1;j<=subspaceSize;j++){
                        minCoeffs[i-1]=minCoeffs[i-1]-H(i,j)*minCoeffs[j-1];
               };

               minCoeffs[i-1]=minCoeffs[i-1]/H(i,i);
        };

};

void GMRES::updateSoln(Grid<wp>&var){
        
        for(int i=0;i<subspaceSize;i++){
                krylovVectors[i]*=minCoeffs[i];
                var+=krylovVectors[i];
        };

};

void GMRES::doIteration(Grid<wp>& var){

        computeResidual(var);

        beta=residualVec.norm2();
        residualVec/=beta;
        krylovVectors[0]=residualVec;

        for(int i=1;i<=subspaceSize;i++){
                ATimesVec(krylovVectors[i-1],krylovVectors[i]);
                GrammSchmidt(i);
/*                if(H(i+1,i)==0.0){
                        subspaceSize=i;
                        break;
                };*/
                krylovVectors[i]/=H(i+1,i);
        };

        givensRotation();
        solveForCoeffs();
        updateSoln(var);

};

void GMRES::solve(Grid<wp>& var){
        currIter=0;

        logger.log("GMRES iteration start: "+std::to_string(residualCalc(var))+"\n",1);

        doIteration(var);

        logger.log("GMRES iteration finished: "+std::to_string(residualCalc(var))+"\n",1);
};

wp GMRES::residualCalc(Grid<wp>& var){
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

bool GMRES::isConverged(Grid<wp>& var){
        if(currIter>maxIters || residualCalc(var)<tol) return true;
        return false;
};
