#include"scheme.h"

BCType StringToBCType(std::string str){
        if(str=="periodic") return periodic;

        if(str=="dirichlet") return dirichlet;

        if(str=="neumann") return neumann;

        if(str=="extrapolate") return extrapolate;

        logger.log("Invalid BC type",1);
        return invalidBCType; 
};

schemeKey StringToschemeKey(std::string str){
        if(str=="timeStepping") return timeStepping;
        
        if(str=="timeSteppingLocality") return timeSteppingLocality;

        if(str=="derivativeX") return derivativeX;

        if(str=="derivativeY") return derivativeY;

        if(str=="derivativeXX") return derivativeXX;

        if(str=="derivativeYY") return derivativeYY;

        if(str=="derivativeXY") return derivativeXY;

        if(str=="faceReconstruction") return faceReconstruction;

        if(str=="stateReconstruction") return stateReconstruction;

        if(str=="fluxCalc") return fluxCalc;

        if(str=="multigrid") return multigrid;

        if(str=="PVCoupling") return PVCoupling;

        if(str=="interfaceConfig") return interfaceConfig;

        logger.log("Invalid scheme key",1);
        return invalidSchemeKey;
};

schemeVal StringToschemeVal(std::string str){
        if(str=="EULER") return EULER;

        if(str=="RK4") return RK4;

        if(str=="LOCAL") return LOCAL;

        if(str=="GLOBAL") return GLOBAL;

        if(str=="LINEARSOLVER") return LINEARSOLVER;

        if(str=="FOU") return FOU;

        if(str=="QUICK") return QUICK;

        if(str=="C4") return C4;

        if(str=="V_CYCLE") return V_CYCLE;

        if(str=="W_CYCLE") return W_CYCLE;

        if(str=="NONE") return NONE;

        if(str=="MUSCL") return MUSCL;

        if(str=="AUSM") return AUSM;

        if(str=="SIMPLE") return SIMPLE;

        if(str=="PISO") return PISO;

        if(str=="IBM") return IBM;

        if(str=="VOF") return VOF;

        if(str=="LEVELSET") return LEVELSET;

        logger.log("Invalid scheme val",1);
        return INVALIDSCHEMEVAL;
};

schemeVal& Scheme::getScheme(schemeKey keyIn){
        return schemeDict[keyIn];
};


void Scheme::setScheme(schemeKey keyIn, schemeVal valIn){

        schemeDict[keyIn]=valIn;
        logger.log("Debug: Scheme setting",1);
};

void Scheme::setBC(std::string boundary, BCType type, wp val){

        BCDict[boundary]={type,val};

        logger.log("Debug: BC setting",1);
};

std::pair<BCType,wp> Scheme::getBC(std::string boundary){
        return BCDict[boundary];

};

bool Scheme::hasScheme(schemeKey key){
        return schemeDict.hasKey(key); 
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

void Scheme::EBD1(Grid<wp>& phi, Grid<wp>& derivative, Mesh& mesh, char dirFlag){

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

void Scheme::ECD2NonUniform(Grid<wp>& phi, Grid<wp>& derivative, Mesh& mesh, char dirFlag){

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

void Scheme::ECD2(Grid<wp>& phi, Grid<wp>& derivative, Mesh& mesh, char dirFlag){

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

void Scheme::EFD1(Grid<wp>& phi, Grid<wp>& derivative, Mesh& mesh, char dirFlag){

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

void Scheme::fillCoeffs(Grid<wp>& phi, char dirFlag){

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

void Scheme::ICD4(Grid<wp>& phi, Grid<wp>& derivative, Mesh& mesh, char dirFlag){

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

Vec2 Scheme::greenGaussCellBased(Grid<wp>&var,Mesh&mesh, int i, int j){

        Vec2 currGrad(0.0,0.0);

        if(i-1>=1){
                currGrad=currGrad-(mesh.normalsI(i,j)*averageValue<wp>(var(i-1,j),var(i,j),mesh,i,j,'x'));
        }
        else{
                currGrad=currGrad-(mesh.normalsI(i,j)*var(i,j));
        };

        if(i+1<=var.intVect[0]){
                currGrad=currGrad+(mesh.normalsI(i+1,j)*averageValue<wp>(var(i,j),var(i+1,j),mesh,i+1,j,'x'));
        }
        else{
                currGrad=currGrad+(mesh.normalsI(i+1,j)*var(i,j));
        };

        if(j-1>=1){
                currGrad=currGrad-(mesh.normalsJ(i,j)*averageValue<wp>(var(i,j-1),var(i,j),mesh,i,j,'y'));
        }
        else{
                currGrad=currGrad-(mesh.normalsJ(i,j)*var(i,j));
        };

        if(j+1<=var.intVect[1]){
                currGrad=currGrad+(mesh.normalsJ(i,j+1)*averageValue<wp>(var(i,j),var(i,j+1),mesh,i,j+1,'y'));
        }
        else{
                currGrad=currGrad+(mesh.normalsJ(i,j+1)*var(i,j));
        };
                        
        currGrad/=mesh.volumes(i,j);

        return currGrad;
};

Vec2 Scheme::adjustedFaceGradient(Grid<wp>&var, Mesh& mesh, int i, int j,char dir){

        Vec2 gradL,gradR,grad;

        gradL=dir=='x'?greenGaussCellBased(var,mesh,i-1,j):greenGaussCellBased(var,mesh,i,j-1);
        gradR=greenGaussCellBased(var,mesh,i,j);

        grad=averageValue<Vec2>(gradL,gradR,mesh,i,j,dir);

        Node ccL,ccR;

        ccL=dir=='x'?mesh.cc(i-1,j):mesh.cc(i,j-1);
        ccR=mesh.cc(i,j);

        Vec2 ccDir(ccR.x-ccL.x,ccR.y-ccL.y);
        wp d=ccDir.norm2();
        ccDir/=d;
        
        wp varL,varR;
        varL=dir=='x'?var(i-1,j):var(i,j-1);
        varR=var(i,j);

        grad=grad+ccDir*((varR-varL)/d-grad.dotProduct(ccDir));

        return grad;
};

std::vector<std::pair<int,wp>>Scheme::FOUReconstruction(wp& massFlux){
        
        std::vector<std::pair<int,wp>>ans;

        if(massFlux>=0.0){
                ans.push_back({-1,1.0});
        }
        else{
                ans.push_back({0,1.0});
        };

        return ans;
};

std::vector<std::pair<int,wp>>Scheme::QUICKReconstruction(wp& massFlux){

        std::vector<std::pair<int,wp>>ans;

        if(massFlux>=0.0){
                ans.push_back({-1,6.0/8.0});
                ans.push_back({0,3.0/8.0});
                ans.push_back({-2,-1.0/8.0});
        }
        else{
                ans.push_back({0,6.0/8.0});
                ans.push_back({-1,3.0/8.0});
                ans.push_back({1,-1.0/8.0});
        };

        return ans;

};

void Scheme::defaultReconstruction(Grid<wp>&var, Grid<wp>&stateL,Grid<wp>&stateR, int i, int j, std::pair<int,int>flags){

        stateL(i,j)=var(i-flags.first,j-flags.second);
        stateR(i,j)=var(i,j);

};


//Flux calc schemes
void Scheme::ausmFlux(dictionary<std::string,wp>& params, std::pair<bool, bool>zeroFluxBC, Grid<Vec2>&normals,
dictionary<std::string,Grid<wp>>&flux, dictionary<std::string,Grid<wp>>&stateL, dictionary<std::string,Grid<wp>>&stateR,
std::pair<int,int>flags){

        Vec2 nHat;
        wp S;
        wp contraVelL,contraVelR,cL,cR,cM;
        wp ML,MR;
        wp alphaLPlus,alphaRMinus;
        wp betaL,betaR;
        wp DLPlus,DRMinus;
        wp DL,DR;
        wp cvlPlus,cvlMinus;
        wp h0L,h0R;
        wp gm=params["GM"];
        wp RGas=params["RGAS"];

        Grid<wp>&rhoL=stateL["rho"];
        Grid<wp>&uL=stateL["u"];
        Grid<wp>&vL=stateL["v"];
        Grid<wp>&pL=stateL["p"];
        Grid<wp>&rhoR=stateR["rho"];
        Grid<wp>&uR=stateR["u"];
        Grid<wp>&vR=stateR["v"];
        Grid<wp>&pR=stateR["p"];

        Grid<wp>&rhoFlux=flux["rho"];
        Grid<wp>&uFlux=flux["u"];
        Grid<wp>&vFlux=flux["v"];
        Grid<wp>&pFlux=flux["p"];

        for(int j=1;j<=uFlux.intVect[1];j++){
                for(int i=1;i<=uFlux.intVect[0];i++){

                        S=normals(i,j).norm2();
                        nHat=normals(i,j)/S;

                        cL=sqrt(gm*pL(i,j)/rhoL(i,j));
                        cR=sqrt(gm*pR(i,j)/rhoR(i,j));
                        cM=0.5*(cL+cR);

                        contraVelL=nHat.dotProduct(Vec2(uL(i,j),vL(i,j)));
                        contraVelR=nHat.dotProduct(Vec2(uR(i,j),vR(i,j)));
                        ML=contraVelL/cM;
                        MR=contraVelR/cM;
                        
                        alphaLPlus=0.5*(1.0+sign(ML));
                        alphaRMinus=0.5*(1.0-sign(MR));

                        betaL=-std::max((wp)0.0,1.0-floor(fabs(ML)));
                        betaR=-std::max((wp)0.0,1.0-floor(fabs(MR)));

                        DLPlus=0.25*(1.0+ML)*(1.0+ML)*(2.0-ML);
                        DRMinus=0.25*(MR-1.0)*(MR-1.0)*(2.0+MR);

                        DL=alphaLPlus*(1.0+betaL)-betaL*DLPlus;
                        DR=alphaRMinus*(1.0+betaR)-betaR*DRMinus;

                        cvlPlus=alphaLPlus*(1.0+betaL)*ML-0.25*betaL*(1.0+ML)*(1.0+ML);
                        cvlMinus=alphaRMinus*(1.0+betaR)*MR+0.25*betaR*(MR-1.0)*(MR-1.0);

                        h0L=(gm/(gm-1.0))*pL(i,j)/rhoL(i,j);
                        h0L=h0L+0.5*(pow(uL(i,j),2.0)+pow(vL(i,j),2.0));

                        h0R=(gm/(gm-1.0))*pR(i,j)/rhoR(i,j);
                        h0R=h0R+0.5*(pow(uR(i,j),2.0)+pow(vR(i,j),2.0));

                        rhoFlux(i,j)=0.0;
                        uFlux(i,j)=0.0;
                        vFlux(i,j)=0.0;
                        pFlux(i,j)=0.0;

                        if(flags.first==1){
                                if(i==1 && zeroFluxBC.first) goto pf;
                                if(i==uFlux.intVect[0] && zeroFluxBC.second) goto pf;

                        }
                        else if(flags.second==1){
                                if(j==1 && zeroFluxBC.first) goto pf;
                                if(j==uFlux.intVect[1] && zeroFluxBC.second) goto pf;
                        };

                        rhoFlux(i,j)+=S*rhoL(i,j)*cvlPlus*cL;
                        uFlux(i,j)+=S*rhoL(i,j)*cvlPlus*cL*uL(i,j);
                        vFlux(i,j)+=S*rhoL(i,j)*cvlPlus*cL*vL(i,j);
                        pFlux(i,j)+=S*rhoL(i,j)*cvlPlus*cL*h0L;

                        rhoFlux(i,j)+=S*rhoR(i,j)*cvlMinus*cR;
                        uFlux(i,j)+=S*rhoR(i,j)*cvlMinus*cR*uR(i,j);
                        vFlux(i,j)+=S*rhoR(i,j)*cvlMinus*cR*vR(i,j);
                        pFlux(i,j)+=S*rhoR(i,j)*cvlMinus*cR*h0R;

                        pf:
                        uFlux(i,j)=uFlux(i,j)+DL*pL(i,j)*normals(i,j).en.x;
                        vFlux(i,j)=vFlux(i,j)+DR*pR(i,j)*normals(i,j).en.y;

                };

        };


};

