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

        return;
};

void Scheme::setBC(int ind, BCType type, wp val){
        BC[ind]=type;
        BCVal[ind]=val;

};

std::pair<BCType,wp> Scheme::getBC(int ind){
        return {BC[ind],BCVal[ind]};

};

void Scheme::TDMA(int n){

       for(int i=1;i<n;i++){
                RHS[i]=RHS[i]-a[i]*RHS[i-1]/b[i-1];
                b[i]=b[i]-a[i]*c[i-1]/b[i-1];
       };

       x[n-1]=RHS[n-1]/b[n-1];

       for(int i=n-2;i>=0;i--){
                x[i]=(RHS[i]-c[i]*x[i+1])/b[i];
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

