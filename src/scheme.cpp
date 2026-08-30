#include"scheme.h"

schemeVal Scheme::getScheme(schemeKey keyIn){
        schemeVal ans=invalidScheme;

        if(keys.size()==0) return ans;

        int counter=0;

        for(auto key:keys){
                if(key==keyIn){
                        ans=values[counter];
                        return ans;
                };

                counter=counter+1;
        };

        return ans;
};

void Scheme::setScheme(schemeKey keyIn, schemeVal valIn){

        int ind=-1;

        int counter=0;


        if(keys.size()==0) goto skip;


        for(auto key:keys){
                if(key==keyIn){
                        ind=counter;
                        goto skip;
                };

                counter=counter+1;
        };

skip:     if(ind==-1){
                keys.push_back(keyIn);
                values.push_back(valIn);
        }
        else{
                values[ind]=valIn;
        };


        logger.log("Debug: Scheme setting",1);
        return;
};

void Scheme::setBC(int ind, BCType type, wp val){

        BC[ind]=type;
        BCVal[ind]=val;


        logger.log("Debug: BC setting",1);
};

std::pair<BCType,wp> Scheme::getBC(int ind){
        return {BC[ind],BCVal[ind]};

};

void Scheme::TDMA(int n, wp* xIn){

       for(int i=1;i<n;i++){
                RHS[i]=RHS[i]-a[i]*RHS[i-1]/b[i-1];
                b[i]=b[i]-a[i]*c[i-1]/b[i-1];
       };

       xIn[n-1]=RHS[n-1]/b[n-1];

       for(int i=n-2;i>=0;i--){
                xIn[i]=(RHS[i]-c[i]*xIn[i+1])/b[i];
       };

};

void Scheme::EBD1(Grid& phi, Grid& derivative, Mesh& mesh, char dirFlag){

       std::array<int,6> gridDims=phi.size();
       
       int imx=gridDims[0];
       int jmx=gridDims[1];

       wp delta;

       if(dirFlag=='x'){
               for(int j=1;j<=jmx;j++){
                  for(int i=1;i<=imx;i++){
                       delta=mesh(i,j).x-mesh(i-1,j).x;
                       derivative(i,j)=(phi(i,j)-phi(i-1,j))/delta;
                  };
               };
       }
       else if(dirFlag=='y'){
               for(int j=1;j<=jmx;j++){
                  for(int i=1;i<=imx;i++){
                       delta=mesh(i,j).y-mesh(i,j-1).y;
                       derivative(i,j)=(phi(i,j)-phi(i,j-1))/delta;
                  };
               };
       };

};

void Scheme::ECD2NonUniform(Grid& phi, Grid& derivative, Mesh& mesh, char dirFlag){

       std::array<int,6> gridDims=phi.size();
       
       int imx=gridDims[0];
       int jmx=gridDims[1];

       wp delta1,delta2,alpha;

       if(dirFlag=='x'){
               for(int j=1;j<=jmx;j++){
                  for(int i=1;i<=imx;i++){
                       delta1=mesh(i,j).x-mesh(i-1,j).x;
                       delta2=mesh(i+1,j).x-mesh(i,j).x;
                       alpha=delta2/delta1;
                       derivative(i,j)=(-phi(i,j)*(1.0-alpha*alpha)+phi(i+1,j)-alpha*alpha*phi(i-1,j))/
                                                                                (alpha*(1.0+alpha)*delta1);
                  };
               };
       }
       else if(dirFlag=='y'){
       
       };

};

void Scheme::ECD2(Grid& phi, Grid& derivative, Mesh& mesh, char dirFlag){

       std::array<int,6> gridDims=phi.size();
       
       int imx=gridDims[0];
       int jmx=gridDims[1];

       wp delta;

       if(dirFlag=='x'){
               for(int j=1;j<=jmx;j++){
                  for(int i=1;i<=imx;i++){
                       delta=mesh(i+1,j).x-mesh(i-1,j).x;
                       derivative(i,j)=(phi(i+1,j)-phi(i-1,j))/delta;
                  };
               };
       }
       else if(dirFlag=='y'){
               for(int j=1;j<=jmx;j++){
                  for(int i=1;i<=imx;i++){
                       delta=mesh(i,j+1).y-mesh(i,j-1).y;
                       derivative(i,j)=(phi(i,j+1)-phi(i,j-1))/delta;
                  };
               };
       };

};

void Scheme::EFD1(Grid& phi, Grid& derivative, Mesh& mesh, char dirFlag){

       std::array<int,6> gridDims=phi.size();
       
       int imx=gridDims[0];
       int jmx=gridDims[1];

       wp delta;

       if(dirFlag=='x'){
               for(int j=1;j<=jmx;j++){
                  for(int i=1;i<=imx;i++){
                       delta=mesh(i+1,j).x-mesh(i,j).x;
                       derivative(i,j)=(phi(i+1,j)-phi(i,j))/delta;
                  };
               };
       }
       else if(dirFlag=='y'){
               for(int j=1;j<=jmx;j++){
                  for(int i=1;i<=imx;i++){
                       delta=mesh(i,j+1).y-mesh(i,j).y;
                       derivative(i,j)=(phi(i,j+1)-phi(i,j))/delta;
                  };
               };
       };

};

void Scheme::fillCoeffs(Grid& phi, char dirFlag){

       std::array<int,6> gridDims=phi.size();
       
       int imx=gridDims[0];
       int jmx=gridDims[1];

       if(dirFlag=='x'){
                  
                  for(int n=0;n<imx;n++){
                           a[n]=0.25;
                           b[n]=1.0;
                           c[n]=0.25;
                  };

                  a[0]=0.0;
                  c[imx-1]=0.0;

       }
       else if(dirFlag=='y'){
                  
                  for(int n=0;n<jmx;n++){
                           a[n]=0.25;
                           b[n]=1.0;
                           c[n]=0.25;
                  };

                  a[0]=0.0;
                  c[jmx-1]=0.0;

       };

};

void Scheme::ICD4(Grid& phi, Grid& derivative, Mesh& mesh, char dirFlag){

       std::array<int,6> gridDims=phi.size();
       
       int imx=gridDims[0];
       int jmx=gridDims[1];

       wp delta=dirFlag=='x'?mesh(2,1).x-mesh(1,1).x:mesh(1,2).y-mesh(1,1).y; //TODO: Make this general
                                                                              //i need not be x dir

       int sysDims=dirFlag=='x'?imx:jmx;

       if(a==nullptr){
                a=new wp[std::max(imx,jmx)];
                b=new wp[std::max(imx,jmx)];
                c=new wp[std::max(imx,jmx)];
                RHS=new wp[std::max(imx,jmx)];
                x=new wp[std::max(imx,jmx)];
                y=new wp[std::max(imx,jmx)];
       };


       wp multiplicand,dividend;
       
       if(dirFlag=='x'){
               for(int j=1;j<=jmx;j++){

                       fillCoeffs(phi,dirFlag);
                       
                       for(int n=0;n<sysDims;n++){
                                RHS[n] = 1.5*(phi(n+2,j)-phi(n,j))/(2.0*delta);
                       };

                       b[0]=2.0;
                       b[sysDims-1]=1.0+((0.25)/(1.0))*0.25;

                       TDMA(imx, x);
                       
                       fillCoeffs(phi,dirFlag);
                       
                       for(int n=0;n<sysDims;n++){

                                RHS[n] = 0.0;
                       };

                       RHS[0]=-1.0;
                       RHS[sysDims-1]=0.25;
                       b[0]=2.0;
                       b[sysDims-1]=1.0+((0.25)/(1.0))*0.25;

                       TDMA(imx, y);

                       multiplicand=1.0*x[0]-((0.25)/(1.0))*x[sysDims-1];
                       dividend=1.0+1.0*y[0]-((0.25)/(1.0))*y[sysDims-1];

                       for(int n=0;n<sysDims;n++){
                                derivative(n+1,j) = x[n]-y[n]*multiplicand/dividend;
                       };

               };
        }
       else{

       };

};

