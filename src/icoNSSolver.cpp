#include"icoNSSolver.h"

icoNSSolver::icoNSSolver(){
        std::ofstream outfile("res.dat");
        outfile<<""<<std::endl;
        outfile.close();

        std::ifstream infile(meshFile);
        int imx, jmx;

        infile>>imx>>jmx;

        mesh = new Mesh(Box(imx,jmx,1,1,1,1),{NODE,NODE});

        for(int j=1;j<=jmx;j++){
                for(int i=1;i<=imx;i++){
                        infile>>(*mesh)(i,j).x>>(*mesh)(i,j).y;
                };
        };

        mesh->setGhostNodes();
        mesh->calcAreas();
        mesh->calcVolumes();

        infile.close();
      
        logger.log("Allocating memory",1);
        (*this).setVar("u",{CELL,CELL});
        (*this).setVar("v",{CELL,CELL});
        (*this).setVar("p",{CELL,CELL});
        (*this).setVar("gradPX",{CELL,CELL});
        (*this).setVar("gradPY",{CELL,CELL});
        (*this).setVar("p'",{CELL,CELL});
        RHS["u"]=vars["u"];
        RHS["v"]=vars["v"];
        RHS["p'"]=vars["p'"];
        coeffs["u"]=Grid<boundMatRow>(Box(imx,jmx,1,1,1,1),{CELL,CELL});
        coeffs["v"]=Grid<boundMatRow>(Box(imx,jmx,1,1,1,1),{CELL,CELL});
        coeffs["p'"]=Grid<boundMatRow>(Box(imx,jmx,1,1,1,1),{CELL,CELL});
        IFaceMDot=Grid<wp>(Box(imx,jmx,1,1,1,1),{NODE,CELL});
        JFaceMDot=Grid<wp>(Box(imx,jmx,1,1,1,1),{CELL,NODE});

        std::string data,data1,data2;
        wp dataNum;
        std::ifstream infile1(settingsFile);
        
        logger.log("Reading input files",1);
        while(infile1>>data){
                if(data.front()=='#') continue;
                
                if(data=="FLOWPROPERTIES"){
                                logger.log("Reading flow properties",1);
                                while(infile1>>data){
                                        if(data.front()=='#') break;
                                        infile1>>dataNum;
                                        params[data]=dataNum;
                                };
                                continue;
                }
                else if(data=="SCHEMES"){
                                logger.log("Reading schemes",1);
                                while(infile1>>data){
                                        if(data.front()=='#') break;
                                        infile1>>data1;
                                        setScheme(StringToschemeKey(data),StringToschemeVal(data1));
                                        if(data1=="LINEARSOLVER"){ 
                                                infile1>>data2;
                                                switch(stringToSolverType(data2)){
                                                        case(JACOBI_SOLVER):
                                                                ls=new Jacobi(100,params["TOL"],1.0,vars["p'"],imx-1,jmx-1); 
                                                                break;
                                                        case(GAUSS_SEIDEL_SOLVER):
                                                                ls=new GaussSeidel(100,params["TOL"],1.0,vars["p'"],imx-1,jmx-1); 
                                                                break;
                                                        case(BICONJUGATE_GRADIENT_SOLVER):
                                                                ls=new biCGStab(100,params["TOL"],1.0,vars["p'"],imx-1,jmx-1);
                                                                break;
                                                        default:
                                                                logger.log("Invalid linear solver name",1);
                                                                break;
                                                };
                                        };
                                };
                                continue;
               }
               else if(data=="BOUNDARYCONDITIONS"){
                                logger.log("Reading BC",1);
                                while(infile1>>data){
                                        if(data.front()=='#') break;
                                        infile1>>data1;
                                        if(data1=="NOSLIPWALL"){
                                                globalBCDict[data]=data1;
                                                infile1>>dataNum;
                                                scheme.BCDict["u"+data]={dirichlet,dataNum};
                                                infile1>>dataNum;
                                                scheme.BCDict["v"+data]={dirichlet,dataNum};
                                                infile1>>dataNum;
                                                scheme.BCDict["p"+data]={neumann,0.0};
                                        }
                                        else if(data1=="INLET"){
                                                globalBCDict[data]=data1;
                                                infile1>>dataNum;
                                                scheme.BCDict["u"+data]={dirichlet,dataNum};
                                                infile1>>dataNum;
                                                scheme.BCDict["v"+data]={dirichlet,dataNum};
                                                infile1>>dataNum;
                                                scheme.BCDict["p"+data]={dirichlet,dataNum};
                                        }
                                        else if(data1=="OUTLET"){
                                                globalBCDict[data]=data1;
                                                infile1>>dataNum;
                                                scheme.BCDict["u"+data]={neumann,0.0};
                                                infile1>>dataNum;
                                                scheme.BCDict["v"+data]={neumann,0.0};
                                                infile1>>dataNum;
                                                scheme.BCDict["p"+data]={neumann,0.0};
                                        }

                                        else{
                                                logger.log("Invalid BC type",1);
                                        };
                                };
                                continue;
               }
               else{
                                std::cout<<"Warning,invalid input file line\n";
                                break;
               };
        };
        infile1.close();


};

void icoNSSolver::initialCondition(){

        std::ifstream infile(initFile);
        std::string temp;
        infile>>temp;
        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        infile>>vars["u"](i,j);
                };
        };

        infile>>temp;
        for(int j=1;j<=vars["v"].intVect[1];j++){
                for(int i=1;i<=vars["v"].intVect[0];i++){
                        infile>>vars["v"](i,j);
                };
        };

        infile>>temp;
        for(int j=1;j<=vars["p"].intVect[1];j++){
                for(int i=1;i<=vars["p"].intVect[0];i++){
                        infile>>vars["p"](i,j);
                };
        };

        infile.close();

        for(int j=1;j<=IFaceMDot.intVect[1];j++){
                for(int i=1;i<=IFaceMDot.intVect[0];i++){
                        IFaceMDot(i,j)=params["RHO"]*mesh->normalsI(i,j).dotProduct(Vec2(0.5*(vars["u"](i-1,j)+vars["u"](i,j)),
                                                                                   0.5*(vars["v"](i-1,j)+vars["v"](i,j))));
                };
        };

        for(int j=1;j<=JFaceMDot.intVect[1];j++){
                for(int i=1;i<=JFaceMDot.intVect[0];i++){
                        JFaceMDot(i,j)=params["RHO"]*mesh->normalsJ(i,j).dotProduct(Vec2(0.5*(vars["u"](i,j-1)+vars["u"](i,j)),
                                                                                   0.5*(vars["v"](i,j-1)+vars["v"](i,j))));
                };
        };

        if(scheme.getScheme(PVCoupling)==PISO){
                tempVar=vars["u"];
        };
};

void icoNSSolver::applyBC(){

        dictionary<std::string,std::pair<int,int>>lo;
        dictionary<std::string,std::pair<int,int>>hi;
        dictionary<std::string,std::pair<std::pair<int,int>,std::pair<int,int>>>inc;
        dictionary<std::string,wp>sig;
     
        lo["BCW"]={1,1};
        hi["BCW"]={1,mesh->box.jmx-1};
        inc["BCW"]={{0,0},{0,1}};
        sig["BCW"]=-1.0;

        lo["BCE"]={mesh->box.imx-1,1};
        hi["BCE"]={mesh->box.imx-1,mesh->box.jmx-1};
        inc["BCE"]={{1,0},{1,1}};
        sig["BCE"]=1.0;

        lo["BCS"]={1,1};
        hi["BCS"]={mesh->box.imx-1,1};
        inc["BCS"]={{0,0},{1,0}};
        sig["BCS"]=-1.0;

        lo["BCN"]={1,mesh->box.jmx-1};
        hi["BCN"]={mesh->box.imx-1,mesh->box.jmx-1};
        inc["BCN"]={{0,1},{1,1}};
        sig["BCN"]=1.0;

        wp mu=params["MU"];
        
        for(auto key:globalBCDict.getKeys()){
                if(globalBCDict[key]=="NOSLIPWALL"){
                        wp uWall,vWall;
                        wp nx,ny,S,d;
                        Vec2 l;

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){
                                        uWall=scheme.BCDict["u"+key].second;
                                        vWall=scheme.BCDict["v"+key].second;

                                        if(key=="BCE" || key=="BCW"){
                                                nx=mesh->normalsI(i+inc[key].first.first,j+inc[key].first.second).en.x;
                                                ny=mesh->normalsI(i+inc[key].first.first,j+inc[key].first.second).en.y;
                                                IFaceMDot(i+inc[key].first.first,j+inc[key].first.second)=0.0;
                                        }
                                        else if(key=="BCS" || key=="BCN"){
                                                nx=mesh->normalsJ(i+inc[key].first.first,j+inc[key].first.second).en.x;
                                                ny=mesh->normalsJ(i+inc[key].first.first,j+inc[key].first.second).en.y;
                                                JFaceMDot(i+inc[key].first.first,j+inc[key].first.second)=0.0;
                                        };

                                        S=sqrt(nx*nx+ny*ny);
                                        nx=nx/S;
                                        ny=ny/S;

                                        l=Vec2((*mesh)(i+inc[key].first.first,j+inc[key].first.second),
                                               (*mesh)(i+inc[key].second.first,j+inc[key].second.second));

                                        d=Mesh::pDist(l,mesh->cc(i,j),Vec2(nx,ny));
                                       
                                        //Momentum eqn
                                        changeMatElement(i,j,i,j,mu*(1.0-nx*nx)*S/d,coeffs["u"]);
                                        RHS["u"](i,j)+=(mu*(uWall*(1.0-nx*nx)+(vars["v"](i,j)-vWall)*nx*ny)*S/d);

                                        changeMatElement(i,j,i,j,mu*(1.0-ny*ny)*S/d,coeffs["v"]);
                                        RHS["v"](i,j)+=(mu*(vWall*(1.0-ny*ny)+(vars["u"](i,j)-uWall)*nx*ny)*S/d);


                                };
                        };

                };

                if(globalBCDict[key]=="INLET"){

                        wp rho=params["RHO"];
                        wp mDot;
                        wp uFree,vFree;
                        wp nx,ny,S;
                        Vec2 l;
                        wp d;
                        wp mu=params["MU"];

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){

                                         //Convective flux
                                         uFree=scheme.BCDict["u"+key].second;
                                         vFree=scheme.BCDict["v"+key].second;
                                        
                                        if(key=="BCW" || key=="BCE"){
                                              mDot=rho*sig[key]*(mesh->normalsI(i+inc[key].first.first,j+inc[key].first.second).
                                                                         dotProduct(Vec2(uFree,vFree)));
                                              IFaceMDot(i+inc[key].first.first,j+inc[key].first.second)=sig[key]*mDot;
                                        }
                                        else if(key=="BCS" || key=="BCN"){
                                              mDot=rho*sig[key]*(mesh->normalsJ(i+inc[key].first.first,j+inc[key].first.second).
                                                                         dotProduct(Vec2(uFree,vFree)));
                                              JFaceMDot(i+inc[key].first.first,j+inc[key].first.second)=sig[key]*mDot;
                                        };

                                        for(auto key1:{"u","v"}){
                                                 RHS[key1](i,j)-=mDot*scheme.BCDict[key1+key].second;
                                        };

                                        //Diffusive flux

                                        if(key=="BCE" || key=="BCW"){
                                                nx=mesh->normalsI(i+inc[key].first.first,j+inc[key].first.second).en.x;
                                                ny=mesh->normalsI(i+inc[key].first.first,j+inc[key].first.second).en.y;
                                        }
                                        else if(key=="BCS" || key=="BCN"){
                                                nx=mesh->normalsJ(i+inc[key].first.first,j+inc[key].first.second).en.x;
                                                ny=mesh->normalsJ(i+inc[key].first.first,j+inc[key].first.second).en.y;
                                        };

                                        S=sqrt(nx*nx+ny*ny);
                                        nx=nx/S;
                                        ny=ny/S;

                                        l=Vec2((*mesh)(i+inc[key].first.first,j+inc[key].first.second),
                                               (*mesh)(i+inc[key].second.first,j+inc[key].second.second));

                                        d=Mesh::pDist(l,mesh->cc(i,j),Vec2(nx,ny));

                                        changeMatElement(i,j,i,j,mu*S/d,coeffs["u"]);
                                        RHS["u"](i,j)+=mu*uFree*S/d;
   
                                        changeMatElement(i,j,i,j,mu*S/d,coeffs["v"]);
                                        RHS["v"](i,j)+=mu*vFree*S/d;
                                       
                                };
                        };

                };

                if(globalBCDict[key]=="OUTLET"){

                        wp rho=params["RHO"];
                        wp mDot;
                        wp uFree,vFree,pFree;
                        wp Du,Dv,Dc;
                        wp S,nx,ny;
                        Vec2 l;
                        wp d;

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){
                                       
                                        //Convective flux
                                        uFree=vars["u"](i,j);
                                        vFree=vars["v"](i,j);
                                        pFree=scheme.BCDict["p"+key].second;

                                        if(key=="BCW" || key=="BCE"){
                                              mDot=rho*sig[key]*(mesh->normalsI(i+inc[key].first.first,j+inc[key].first.second).
                                                                         dotProduct(Vec2(uFree,vFree)));
                                              IFaceMDot(i+inc[key].first.first,j+inc[key].first.second)=sig[key]*mDot;
                                        }
                                        else if(key=="BCS" || key=="BCN"){
                                              mDot=rho*sig[key]*(mesh->normalsJ(i+inc[key].first.first,j+inc[key].first.second).
                                                                         dotProduct(Vec2(uFree,vFree)));
                                              JFaceMDot(i+inc[key].first.first,j+inc[key].first.second)=sig[key]*mDot;
                                        };

                                        for(auto key1:{"u","v"}){
                                                 changeMatElement(i,j,i,j,mDot,coeffs[key1]);
                                        };
                                         
                                        //Diffusive flux
                                        //Terms go to 0 because of extrapolation


                                        //Pressure in momentum equation: Need to take outlet pressure into account when computing gradP
                                        if(key=="BCE" || key=="BCW"){
                                                nx=mesh->normalsI(i+inc[key].first.first,j+inc[key].first.second).en.x;
                                                ny=mesh->normalsI(i+inc[key].first.first,j+inc[key].first.second).en.y;
                                        }
                                        else if(key=="BCS" || key=="BCN"){
                                                nx=mesh->normalsJ(i+inc[key].first.first,j+inc[key].first.second).en.x;
                                                ny=mesh->normalsJ(i+inc[key].first.first,j+inc[key].first.second).en.y;
                                        };

                                        S=sqrt(nx*nx+ny*ny);
                                        nx=nx/S;
                                        ny=ny/S;

                                        l=Vec2((*mesh)(i+inc[key].first.first,j+inc[key].first.second),
                                               (*mesh)(i+inc[key].second.first,j+inc[key].second.second));

                                        d=Mesh::pDist(l,mesh->cc(i,j),Vec2(nx,ny));

                                        vars["gradPX"](i,j)=sig[key]*(pFree-vars["p"](i,j))*nx/d;
                                        vars["gradPY"](i,j)=sig[key]*(pFree-vars["p"](i,j))*ny/d;

                                        //Pressure correction eqn
                                        Du=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["u"]);
                                        Dv=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["v"]);
                                        if(key=="BCW" || key=="BCE"){
                                                S=mesh->normalsI(i+inc[key].first.first,j+inc[key].first.second).norm2();
                                                Dc=mesh->normalsI(i+inc[key].first.first,j+inc[key].first.second).
                                                   dotProduct(Vec2(Du,Dv))/S;
                                        }
                                        else if(key=="BCS" || key=="BCN"){
                                                S=mesh->normalsJ(i+inc[key].first.first,j+inc[key].first.second).norm2();
                                                Dc=mesh->normalsJ(i+inc[key].first.first,j+inc[key].first.second).
                                                   dotProduct(Vec2(Du,Dv))/S;
                                        };

                                        changeMatElement(i,j,i,j,-rho*Dc,coeffs["p'"]); //TODO: Fix
                                        //TODO: Will outlet get p'
                                };
                        };

                };


        };


        //Stabilize FOU by subtracting out net mass flux

        wp mDot;

        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        mDot=IFaceMDot(i+1,j)-IFaceMDot(i,j)+JFaceMDot(i,j+1)-JFaceMDot(i,j);
                        changeMatElement(i,j,i,j,-mDot,coeffs["u"]);
                        changeMatElement(i,j,i,j,-mDot,coeffs["v"]);
                };
        };

};

void icoNSSolver::QDot(){

};

void icoNSSolver::lSolve(){

        if(scheme.getScheme(PVCoupling)==SIMPLE){
                SIMPLEDriver();
        }
        else if(scheme.getScheme(PVCoupling)==PISO){
                PISODriver();
        }
        else{
                logger.log("Invalid PV coupling algorithm specified; solver terminating",1);
                std::exit(-1);
        };

        getResNorm();
};

void icoNSSolver::computeTimeStep(Grid<wp>&dt){

};

void icoNSSolver::updateVars(Grid<wp>&dt, wp storeFactor){

};

void icoNSSolver::updateVars(Grid<wp>&dt, wp storeFactor, std::vector<Grid<wp>>& uStore){

};

void icoNSSolver::updateVars(Grid<wp>&dt, wp storeFactor, std::vector<Grid<wp>>& uStore, std::vector<Grid<wp>>& rStore){

};

wp icoNSSolver::getResNorm(){

        wp res=0.0;
        wp tempRes;

        std::ofstream outfile("res.dat",std::ios::app);

        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                       tempRes=RHS["u"](i,j);
                       for(auto c:coeffs["u"](i,j)){
                                tempRes=tempRes-c.second*vars["u"](c.first.first,c.first.second);
                       };
                        
                       res+=pow(tempRes,2.0);
                };
        };

        resXMom=res;
       
        std::cout<<"X-mom res: "<<std::scientific<<res<<std::endl;

        res=0.0;
        for(int j=1;j<=vars["v"].intVect[1];j++){
                for(int i=1;i<=vars["v"].intVect[0];i++){
                       tempRes=RHS["v"](i,j);
                       for(auto c:coeffs["v"](i,j)){
                                tempRes=tempRes-c.second*vars["v"](c.first.first,c.first.second);
                       };
                        
                       res+=pow(tempRes,2.0);
                };
        };

        resYMom=res;
        std::cout<<"Y-mom res: "<<std::scientific<<res<<std::endl;

        res=0.0;
        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                       tempRes=IFaceMDot(i+1,j)-IFaceMDot(i,j)+JFaceMDot(i,j+1)-JFaceMDot(i,j);
                       res+=pow(tempRes,2.0);
                };
        };

        resMass=res;
        std::cout<<"Mass balance: "<<std::scientific<<res<<std::endl;

        outfile<<resXMom<<"\t"<<resYMom<<"\t"<<resMass<<std::endl;
        outfile.close();
};

void icoNSSolver::getResidual(){

};

bool icoNSSolver::isConverged(){
        if(resXMom<params["TOL"] && resYMom<params["TOL"] && resMass<params["TOL"]) return true;
        return false;
};

void icoNSSolver::SIMPLEDriver(){
        //Reset coefficients
        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        RHS["u"](i,j)=0.0;
                        for(auto& l:coeffs["u"](i,j)){
                                l.second=0.0;
                        };

                        RHS["v"](i,j)=0.0;
                        for(auto& l:coeffs["v"](i,j)){
                                l.second=0.0;
                        };

                        RHS["p'"](i,j)=0.0;
                        for(auto& l:coeffs["p'"](i,j)){
                                l.second=0.0;
                        };

                };
        };


        logger.log("Cleared coeffs",1);

        //Compute coefficients
        computeDiffusiveFlux("u");
        computeDiffusiveFlux("v");

        computeConvectiveFlux();
        
        computeStoreGrad(vars["p"]);

        //Apply BC
        applyBC();

        //Under-relax
        wp ap;
        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        ap=findMatElement(i,j,i,j,coeffs["u"]);
                        changeMatElement(i,j,i,j,ap*(1.0-alpha)/alpha,coeffs["u"]);
                        RHS["u"](i,j)+=(1.0-alpha)*ap*vars["u"](i,j)/alpha;
                };
        };

        for(int j=1;j<=vars["v"].intVect[1];j++){
                for(int i=1;i<=vars["v"].intVect[0];i++){
                        ap=findMatElement(i,j,i,j,coeffs["v"]);
                        changeMatElement(i,j,i,j,ap*(1.0-alpha)/alpha,coeffs["v"]);
                        RHS["v"](i,j)+=(1.0-alpha)*ap*vars["v"](i,j)/alpha;
                };
        };

        addPressureSource();

        logger.log("Solving x momentum",1);
        //Solve x momentum
        ls->setMatFunc(coeffs["u"]);
        ls->setRHSFunc(RHS["u"]);
        ls->solve(vars["u"]);

        logger.log("Solving y momentum",1);
        //Solve y momentum
        ls->setMatFunc(coeffs["v"]);
        ls->setRHSFunc(RHS["v"]);
        ls->solve(vars["v"]);

        RhieChow();
        assemblePressureCorrectionEqn();

        logger.log("Solving pressure correction",1);

        //Solve pressure-correction eqn
        vars["p'"].initVal(0.0); 
        ls->setMatFunc(coeffs["p'"]);
        ls->setRHSFunc(RHS["p'"]);
        ls->solve(vars["p'"]);

        correctPressure();
        correctMassFlux();
        correctVelocity();

        std::ofstream outfile1("outu.dat");
        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        outfile1<<mesh->cc(i,j).x<<"\t"<<mesh->cc(i,j).y<<"\t"<<vars["u"](i,j)<<std::endl; 
                 };
                 outfile1<<std::endl;
        };

        outfile1.close();

        std::ofstream outfile2("outv.dat");
        for(int j=1;j<=vars["v"].intVect[1];j++){
                for(int i=1;i<=vars["v"].intVect[0];i++){
                        outfile2<<mesh->cc(i,j).x<<"\t"<<mesh->cc(i,j).y<<"\t"<<vars["v"](i,j)<<std::endl; 
                 };
                 outfile2<<std::endl;
        };

        outfile2.close();

        std::ofstream outfile3("outp.dat");
        for(int j=1;j<=vars["p"].intVect[1];j++){
                for(int i=1;i<=vars["p"].intVect[0];i++){
                        outfile3<<mesh->cc(i,j).x<<"\t"<<mesh->cc(i,j).y<<"\t"<<vars["p"](i,j)<<std::endl; 
                 };
                 outfile3<<std::endl;
        };

        outfile3.close();

};      

void icoNSSolver::PISODriver(){
        
        //First, call SIMPLE
        SIMPLEDriver();

        //Reset coefficients
        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        RHS["u"](i,j)=0.0;
                        for(auto& l:coeffs["u"](i,j)){
                                l.second=0.0;
                        };

                        RHS["v"](i,j)=0.0;
                        for(auto& l:coeffs["v"](i,j)){
                                l.second=0.0;
                        };

                        RHS["p'"](i,j)=0.0;
                        for(auto& l:coeffs["p'"](i,j)){
                                l.second=0.0;
                        };

                };
        };

        //Compute coefficients
        computeDiffusiveFlux("u");
        computeDiffusiveFlux("v");

        computeConvectiveFlux();
        
        computeStoreGrad(vars["p"]);

        //Apply BC
        applyBC();

        //Under-relax
        wp ap;
        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        ap=findMatElement(i,j,i,j,coeffs["u"]);
                        changeMatElement(i,j,i,j,ap*(1.0-alpha)/alpha,coeffs["u"]);
                        RHS["u"](i,j)+=(1.0-alpha)*ap*vars["u"](i,j)/alpha;
                };
        };

        for(int j=1;j<=vars["v"].intVect[1];j++){
                for(int i=1;i<=vars["v"].intVect[0];i++){
                        ap=findMatElement(i,j,i,j,coeffs["v"]);
                        changeMatElement(i,j,i,j,ap*(1.0-alpha)/alpha,coeffs["v"]);
                        RHS["v"](i,j)+=(1.0-alpha)*ap*vars["v"](i,j)/alpha;
                };
        };

        addPressureSource();

        //Solve momentum equation(s) explicitly
        wp diag;
        
        tempVar=vars["u"];

        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        vars["u"](i,j)=RHS["u"](i,j);
                        for(auto c:coeffs["u"](i,j)){
                                if(c.first.first==i && c.first.second==j){
                                        diag=c.second;
                                        continue;
                                };

                                vars["u"](i,j)=vars["u"](i,j)-c.second*tempVar(c.first.first,c.first.second);
                        };

                        vars["u"](i,j)/=diag;

                };
        };
       
        tempVar=vars["v"];

        for(int j=1;j<=vars["v"].intVect[1];j++){
                for(int i=1;i<=vars["v"].intVect[0];i++){
                        vars["v"](i,j)=RHS["v"](i,j);
                        for(auto c:coeffs["v"](i,j)){
                                if(c.first.first==i && c.first.second==j){
                                        diag=c.second;
                                        continue;
                                };

                                vars["v"](i,j)=vars["v"](i,j)-c.second*tempVar(c.first.first,c.first.second);
                        };

                        vars["v"](i,j)/=diag;

                };
        };


        //Rhie-Chow interpolation
        RhieChow();

        //Assemble pressure correction eqn
        assemblePressureCorrectionEqn();

        //Solve pressure-correction eqn
        vars["p'"].initVal(0.0);
        ls->setMatFunc(coeffs["p'"]);
        ls->setRHSFunc(RHS["p'"]);
        ls->solve(vars["p'"]);

        correctPressure();
        correctMassFlux();
        correctVelocity();

};

void icoNSSolver::computeDiffusiveFlux(std::string varName){

        Vec2 grad;
        wp d;
        Node pL,pR;
        Vec2 vCC,vCD;
        wp E,T;
        wp gradCD;

        wp mu=params["MU"];

        for(int j=1;j<=IFaceMDot.intVect[1];j++){
                for(int i=2;i<IFaceMDot.intVect[0];i++){
                        pL=mesh->cc(i-1,j);
                        pR=mesh->cc(i,j);
                        vCC=Vec2(pR.x-pL.x,pR.y-pL.y);
                        d=vCC.norm2();
                        E=mesh->normalsI(i,j).dotProduct(vCC)/d;
                        vCD=mesh->normalsI(i,j)-vCC*(E/d);
                        T=vCD.norm2();
                        grad=scheme.adjustedFaceGradient(vars[varName],*mesh,i,j,'x');
                        gradCD=grad.dotProduct(vCD);

                        //Update cell coefficients of i-1,j
                        changeMatElement(i-1,j,i-1,j,mu*E/d,coeffs[varName]);
                        changeMatElement(i-1,j,i,j,-mu*E/d,coeffs[varName]);
                        RHS[varName](i-1,j)+=gradCD*mu;
                        
                        //Update cell coefficients of i,j
                        changeMatElement(i,j,i,j,mu*E/d,coeffs[varName]);
                        changeMatElement(i,j,i-1,j,-mu*E/d,coeffs[varName]);
                        RHS[varName](i,j)-=gradCD*mu;

                };
        };


        for(int j=2;j<JFaceMDot.intVect[1];j++){
                for(int i=1;i<=JFaceMDot.intVect[0];i++){
                        pL=mesh->cc(i,j-1);
                        pR=mesh->cc(i,j);
                        vCC=Vec2(pR.x-pL.x,pR.y-pL.y);
                        d=vCC.norm2();
                        E=mesh->normalsJ(i,j).dotProduct(vCC)/d;
                        vCD=mesh->normalsJ(i,j)-vCC*(E/d);
                        T=vCD.norm2();
                        grad=scheme.adjustedFaceGradient(vars[varName],*mesh,i,j,'y');
                        gradCD=grad.dotProduct(vCD);

                        //Update cell coefficients of i,j-1
                        changeMatElement(i,j-1,i,j-1,mu*E/d,coeffs[varName]);
                        changeMatElement(i,j-1,i,j,-mu*E/d,coeffs[varName]);
                        RHS[varName](i,j-1)+=gradCD*mu;
                        
                        //Update cell coefficients of i,j
                        changeMatElement(i,j,i,j,mu*E/d,coeffs[varName]);
                        changeMatElement(i,j,i,j-1,-mu*E/d,coeffs[varName]);
                        RHS[varName](i,j)-=gradCD*mu;

                };
        };

};

void icoNSSolver::computeConvectiveFlux(){

        wp rho=params["RHO"];

        std::vector<std::pair<int,wp>>recon;

        wp uRecon,vRecon;

        for(int j=1;j<=IFaceMDot.intVect[1];j++){
                for(int i=2;i<IFaceMDot.intVect[0];i++){

                        recon=scheme.FOUReconstruction(IFaceMDot(i,j));
                        

                        for(auto key:{"u","v"}){
                               for(auto r:recon){ 
                                        changeMatElement(i-1,j,i+r.first,j,IFaceMDot(i,j)*r.second,coeffs[key]);
                                        changeMatElement(i,j,i+r.first,j,-IFaceMDot(i,j)*r.second,coeffs[key]);
                               };

                        };

                        if(scheme.getScheme(faceReconstruction)==QUICK && i>2 && i<IFaceMDot.intVect[0]-1){
                              
                              for(auto key:{"u","v"}){
                                        for(auto r:recon){       
                                                RHS[key](i-1,j)+=vars[key](i+r.first,j)*IFaceMDot(i,j)*r.second;
                                                RHS[key](i,j)-=vars[key](i+r.first,j)*IFaceMDot(i,j)*r.second;
                                        };

                              };

                              recon=scheme.QUICKReconstruction(IFaceMDot(i,j));

                              for(auto key:{"u","v"}){
                                        for(auto r:recon){       
                                                RHS[key](i-1,j)-=vars[key](i+r.first,j)*IFaceMDot(i,j)*r.second;
                                                RHS[key](i,j)+=vars[key](i+r.first,j)*IFaceMDot(i,j)*r.second;
                                        };

                              };

                        };

                };
        };

        for(int j=2;j<JFaceMDot.intVect[1];j++){
                for(int i=1;i<=JFaceMDot.intVect[0];i++){

                        recon=scheme.FOUReconstruction(JFaceMDot(i,j));
                        for(auto key:{"u","v"}){
                               for(auto r:recon){ 
                                        changeMatElement(i,j-1,i,j+r.first,JFaceMDot(i,j)*r.second,coeffs[key]);
                                        changeMatElement(i,j,i,j+r.first,-JFaceMDot(i,j)*r.second,coeffs[key]);
                               };
                        };


                        if(scheme.getScheme(faceReconstruction)==QUICK && j>2 && j<JFaceMDot.intVect[1]-1){
                              
                              for(auto key:{"u","v"}){
                                        for(auto r:recon){       
                                                RHS[key](i,j-1)+=vars[key](i,j+r.first)*JFaceMDot(i,j)*r.second;
                                                RHS[key](i,j)-=vars[key](i,j+r.first)*JFaceMDot(i,j)*r.second;
                                        };

                              };

                              recon=scheme.QUICKReconstruction(JFaceMDot(i,j));
                              
                              for(auto key:{"u","v"}){
                                        for(auto r:recon){       
                                                RHS[key](i,j-1)-=vars[key](i,j+r.first)*JFaceMDot(i,j)*r.second;
                                                RHS[key](i,j)+=vars[key](i,j+r.first)*JFaceMDot(i,j)*r.second;
                                        };

                              };

                        };

                };
        };


};

void icoNSSolver::computeStoreGrad(Grid<wp>&var){

        Vec2 grad;
        wp gradPX,gradPY;
        wp state;

        Vec2 dirX,dirY;

        dirX=Vec2(1.0,0.0);
        dirY=Vec2(0.0,1.0);

        for(int j=1;j<=vars["gradPX"].intVect[1];j++){
                for(int i=1;i<=vars["gradPY"].intVect[0];i++){
                        gradPX=0.0;
                        gradPY=0.0;

                        grad=scheme.greenGaussCellBased(var,*mesh,i,j);

                        vars["gradPX"](i,j)=grad.dotProduct(dirX);
                        vars["gradPY"](i,j)=grad.dotProduct(dirY);

                };
        };

};

void icoNSSolver::addPressureSource(){

        Vec2 gradP;
        wp gradPX,gradPY;
        wp au,av;
        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        gradPX=vars["gradPX"](i,j);
                        gradPY=vars["gradPY"](i,j);
                        
                        RHS["u"](i,j)-=gradPX*mesh->volumes(i,j);
                        RHS["v"](i,j)-=gradPY*mesh->volumes(i,j);

                };
        };

};

void icoNSSolver::RhieChow(){

        wp df,dfL,dfR;
        Vec2 gradCorr;
        wp interpU,interpV,interpVelTimesArea;

        wp rho=params["RHO"];

        for(int j=1;j<=IFaceMDot.intVect[1];j++){
                for(int i=2;i<IFaceMDot.intVect[0];i++){

                       gradCorr=scheme.adjustedFaceGradient(vars["p"],*mesh,i,j,'x')-
                                Vec2(scheme.averageValue<wp>(vars["gradPX"](i-1,j),vars["gradPX"](i,j),*mesh,i,j,'x'),
                                     scheme.averageValue<wp>(vars["gradPY"](i-1,j),vars["gradPY"](i,j),*mesh,i,j,'x'));
                       
                       dfL=mesh->volumes(i-1,j)/findMatElement(i-1,j,i-1,j,coeffs["u"]);
                       dfR=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["u"]);
                       df=scheme.averageValue<wp>(dfL,dfR,*mesh,i,j,'x');
                       interpU=scheme.averageValue<wp>(vars["u"](i-1,j),vars["u"](i,j),*mesh,i,j,'x');
                       interpU=interpU-df*gradCorr.dotProduct(Vec2(1.0,0.0));
                       
                       dfL=mesh->volumes(i-1,j)/findMatElement(i-1,j,i-1,j,coeffs["v"]);
                       dfR=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["v"]);
                       df=scheme.averageValue<wp>(dfL,dfR,*mesh,i,j,'x');
                       interpV=scheme.averageValue<wp>(vars["v"](i-1,j),vars["v"](i,j),*mesh,i,j,'x');
                       interpV=interpV-df*gradCorr.dotProduct(Vec2(0.0,1.0));

                       interpVelTimesArea=mesh->normalsI(i,j).dotProduct(Vec2(interpU,interpV));
                       IFaceMDot(i,j)=interpVelTimesArea*rho;

                };
        };

        for(int j=2;j<JFaceMDot.intVect[1];j++){
                for(int i=1;i<=JFaceMDot.intVect[0];i++){

                       gradCorr=scheme.adjustedFaceGradient(vars["p"],*mesh,i,j,'y')-
                                Vec2(scheme.averageValue<wp>(vars["gradPX"](i,j-1),vars["gradPX"](i,j),*mesh,i,j,'y'),
                                     scheme.averageValue<wp>(vars["gradPY"](i,j-1),vars["gradPY"](i,j),*mesh,i,j,'y'));

                       dfL=mesh->volumes(i,j-1)/findMatElement(i,j-1,i,j-1,coeffs["u"]);
                       dfR=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["u"]);
                       df=scheme.averageValue<wp>(dfL,dfR,*mesh,i,j,'y');
                       interpU=scheme.averageValue<wp>(vars["u"](i,j-1),vars["u"](i,j),*mesh,i,j,'y');
                       interpU=interpU-df*gradCorr.dotProduct(Vec2(1.0,0.0));
                       
                       dfL=mesh->volumes(i,j-1)/findMatElement(i,j-1,i,j-1,coeffs["v"]);
                       dfR=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["v"]);
                       df=scheme.averageValue<wp>(dfL,dfR,*mesh,i,j,'y');
                       interpV=scheme.averageValue<wp>(vars["v"](i,j-1),vars["v"](i,j),*mesh,i,j,'y');
                       interpV=interpV-df*gradCorr.dotProduct(Vec2(0.0,1.0));

                       interpVelTimesArea=mesh->normalsJ(i,j).dotProduct(Vec2(interpU,interpV));
                       JFaceMDot(i,j)=interpVelTimesArea*rho;

                };
        };

};

void icoNSSolver::assemblePressureCorrectionEqn(){
       
        wp Ef,Df,d;
        Vec2 dCC;

        Node ccL,ccR;
        wp DLU,DRU,DLV,DRV;

        wp rho=params["RHO"];

        for(int j=1;j<=IFaceMDot.intVect[1];j++){
                for(int i=2;i<IFaceMDot.intVect[0];i++){
                        ccL=mesh->cc(i-1,j);
                        ccR=mesh->cc(i,j);
                        dCC=Vec2(ccR.x-ccL.x,ccR.y-ccL.y);
                        d=dCC.norm2();
                        DLU=mesh->volumes(i-1,j)/findMatElement(i-1,j,i-1,j,coeffs["u"]);
                        DRU=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["u"]);
                        DLV=mesh->volumes(i-1,j)/findMatElement(i-1,j,i-1,j,coeffs["v"]);
                        DRV=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["v"]);

                        Ef=(dCC.en.x*scheme.averageValue<wp>(DLU,DRU,*mesh,i,j,'x')*mesh->normalsI(i,j).en.x+
                           dCC.en.y*scheme.averageValue<wp>(DLV,DRV,*mesh,i,j,'x')*mesh->normalsI(i,j).en.y)/d;

                        Df=Ef/d;
                        changeMatElement(i-1,j,i-1,j,rho*Df,coeffs["p'"]);
                        changeMatElement(i-1,j,i,j,-rho*Df,coeffs["p'"]);

                        changeMatElement(i,j,i,j,rho*Df,coeffs["p'"]);
                        changeMatElement(i,j,i-1,j,-rho*Df,coeffs["p'"]);

                };
        };

        for(int j=2;j<JFaceMDot.intVect[1];j++){
                for(int i=1;i<=JFaceMDot.intVect[0];i++){

                        ccL=mesh->cc(i,j-1);
                        ccR=mesh->cc(i,j);
                        dCC=Vec2(ccR.x-ccL.x,ccR.y-ccL.y);
                        d=dCC.norm2();
                        DLU=mesh->volumes(i,j-1)/findMatElement(i,j-1,i,j-1,coeffs["u"]);
                        DRU=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["u"]);
                        DLV=mesh->volumes(i,j-1)/findMatElement(i,j-1,i,j-1,coeffs["v"]);
                        DRV=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["v"]);

                        Ef=(dCC.en.x*scheme.averageValue<wp>(DLU,DRU,*mesh,i,j,'y')*mesh->normalsJ(i,j).en.x+
                           dCC.en.y*scheme.averageValue<wp>(DLV,DRV,*mesh,i,j,'y')*mesh->normalsJ(i,j).en.y)/d;

                        Df=Ef/d;
                        changeMatElement(i,j-1,i,j-1,rho*Df,coeffs["p'"]);
                        changeMatElement(i,j-1,i,j,-rho*Df,coeffs["p'"]);

                        changeMatElement(i,j,i,j,rho*Df,coeffs["p'"]);
                        changeMatElement(i,j,i,j-1,-rho*Df,coeffs["p'"]);

                };
        };

        for(int j=1;j<=vars["p'"].intVect[1];j++){
                for(int i=1;i<=vars["p'"].intVect[0];i++){

                        RHS["p'"](i,j)-=(IFaceMDot(i+1,j)-IFaceMDot(i,j)
                                        +JFaceMDot(i,j+1)-JFaceMDot(i,j));
                };
        };

};

void icoNSSolver::correctPressure(){

        wp rel=params["REL"];

        for(int j=1;j<=vars["p"].intVect[1];j++){
                for(int i=1;i<=vars["p"].intVect[0];i++){
                        vars["p"](i,j)+=rel*vars["p'"](i,j);
                };
        };

};

void icoNSSolver::correctMassFlux(){

        wp rho=params["RHO"];
        wp dfL,dfR,df;
        Vec2 gradCorr;
        wp mDotCorrX,mDotCorrY;

        for(int j=1;j<=IFaceMDot.intVect[1];j++){
                for(int i=2;i<IFaceMDot.intVect[0];i++){

                       gradCorr=scheme.adjustedFaceGradient(vars["p'"],*mesh,i,j,'x');

                       dfL=mesh->volumes(i-1,j)/findMatElement(i-1,j,i-1,j,coeffs["u"]);
                       dfR=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["u"]);
                       df=scheme.averageValue<wp>(dfL,dfR,*mesh,i,j,'x');
                       
                       mDotCorrX=rho*df*gradCorr.dotProduct(Vec2(1.0,0.0));

                       dfL=mesh->volumes(i-1,j)/findMatElement(i-1,j,i-1,j,coeffs["v"]);
                       dfR=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["v"]);
                       df=scheme.averageValue<wp>(dfL,dfR,*mesh,i,j,'x');

                       mDotCorrY=rho*df*gradCorr.dotProduct(Vec2(0.0,1.0));

                       IFaceMDot(i,j)-=mesh->normalsI(i,j).dotProduct(Vec2(mDotCorrX,mDotCorrY));
                
                };
        };

        for(int j=2;j<JFaceMDot.intVect[1];j++){
                for(int i=1;i<=JFaceMDot.intVect[0];i++){

                       gradCorr=scheme.adjustedFaceGradient(vars["p'"],*mesh,i,j,'y');

                       dfL=mesh->volumes(i,j-1)/findMatElement(i,j-1,i,j-1,coeffs["u"]);
                       dfR=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["u"]);
                       df=scheme.averageValue<wp>(dfL,dfR,*mesh,i,j,'y');
                       
                       mDotCorrX=rho*df*gradCorr.dotProduct(Vec2(1.0,0.0));

                       dfL=mesh->volumes(i,j-1)/findMatElement(i,j-1,i,j-1,coeffs["v"]);
                       dfR=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["v"]);
                       df=scheme.averageValue<wp>(dfL,dfR,*mesh,i,j,'y');

                       mDotCorrY=rho*df*gradCorr.dotProduct(Vec2(0.0,1.0));

                       JFaceMDot(i,j)-=mesh->normalsJ(i,j).dotProduct(Vec2(mDotCorrX,mDotCorrY));

                };
        };

};

void icoNSSolver::correctVelocity(){

        wp gradPPrimeX,gradPPrimeY;
        wp state;
        Vec2 dirX,dirY;

        dirX=Vec2(1.0,0.0);
        dirY=Vec2(0.0,1.0);

        for(int j=1;j<=vars["p'"].intVect[1];j++){
                for(int i=1;i<=vars["p'"].intVect[0];i++){

                        gradPPrimeX=0.0;
                        gradPPrimeY=0.0;

                        if(i>1){
                                state=scheme.averageValue<wp>(vars["p'"](i-1,j),vars["p'"](i,j),*mesh,i,j,'x');
                        }
                        else{
                                state=vars["p'"](i,j);
                        };

                        gradPPrimeX-=state*mesh->normalsI(i,j).dotProduct(dirX);
                        gradPPrimeY-=state*mesh->normalsI(i,j).dotProduct(dirY);

                        if(i<vars["p'"].intVect[0]){
                                state=scheme.averageValue<wp>(vars["p'"](i,j),vars["p'"](i+1,j),*mesh,i+1,j,'x');
                        }
                        else{
                                state=vars["p'"](i,j);
                        };

                        gradPPrimeX+=state*mesh->normalsI(i+1,j).dotProduct(dirX);
                        gradPPrimeY+=state*mesh->normalsI(i+1,j).dotProduct(dirY);

                        if(j>1){
                                state=scheme.averageValue<wp>(vars["p'"](i,j-1),vars["p'"](i,j),*mesh,i,j,'y');
                        }
                        else{
                                state=vars["p'"](i,j);
                        };


                        gradPPrimeX-=state*mesh->normalsJ(i,j).dotProduct(dirX);
                        gradPPrimeY-=state*mesh->normalsJ(i,j).dotProduct(dirY);

                        if(j<vars["p'"].intVect[1]){
                                state=scheme.averageValue<wp>(vars["p'"](i,j),vars["p'"](i,j+1),*mesh,i,j+1,'y');
                        }
                        else{
                                state=vars["p'"](i,j);
                        };

                        gradPPrimeX+=state*mesh->normalsJ(i,j+1).dotProduct(dirX);
                        gradPPrimeY+=state*mesh->normalsJ(i,j+1).dotProduct(dirY);

                        gradPPrimeX/=mesh->volumes(i,j);
                        gradPPrimeY/=mesh->volumes(i,j);

                        vars["u"](i,j)-=(mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["u"]))*gradPPrimeX;
                        vars["v"](i,j)-=(mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["v"]))*gradPPrimeY;

                };
        };


};

