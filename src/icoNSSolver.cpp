#include"icoNSSolver.h"

icoNSSolver::icoNSSolver(){

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
                                        if(data1=="LINEARSOLVER") infile1>>data2;
                                        setScheme(StringToschemeKey(data),StringToschemeVal(data1));
                                        switch(stringToSolverType(data2)){
                                                case(JACOBI_SOLVER):
                                                        ls=new Jacobi(500,1E-15,0.5,vars["p"],imx-1,jmx-1); 
                                                        break;
                                                case(GAUSS_SEIDEL_SOLVER):
                                                        ls=new GaussSeidel(500,1E-15,0.5,vars["p"],imx-1,jmx-1); 
                                                        break;
                                                default:
                                                        logger.log("Invalid linear solver name",1);
                                                        break;
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

        IFaceMDot.initVal(params["RHO"]*params["UREF"]*mesh->normalsI(1,1).norm2());
        JFaceMDot.initVal(params["RHO"]*params["VREF"]*mesh->normalsJ(1,1).norm2());

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

                                        changeMatElement(i,j,i,j,mu*(1.0-nx*nx)*S/d,coeffs["u"]);
                                        RHS["u"](i,j)+=mu*(uWall*(1.0-nx*nx)+(vars["v"](i,j)-vWall)*nx*ny)*S/d;

                                        changeMatElement(i,j,i,j,mu*(1.0-ny*ny)*S/d,coeffs["v"]);
                                        RHS["v"](i,j)+=mu*(vWall*(1.0-ny*ny)+(vars["u"](i,j)-uWall)*nx*ny)*S/d;

                                };
                        };

                };

                if(globalBCDict[key]=="INLET"){

                        wp rho=params["RHO"];
                        wp mDot;
                        wp uFree,vFree;

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){

                                         uFree=scheme.BCDict["u"+key].second;
                                         vFree=scheme.BCDict["v"+key].second;

                                        if(key=="BCW" || key=="BCE"){
                                              mDot=rho*sig[key]*(mesh->normalsI(i+inc[key].first.first,j+inc[key].first.second).
                                                                         dotProduct(Vec2(uFree,vFree)));
                                        }
                                        else if(key=="BCS" || key=="BCN"){
                                              mDot=rho*sig[key]*(mesh->normalsJ(i+inc[key].first.first,j+inc[key].first.second).
                                                                         dotProduct(Vec2(uFree,vFree)));
                                        };

                                        for(auto key1:vars.getKeys()){
                                                 if(key=="p") continue;
                                                 RHS[key1](i,j)-=mDot*scheme.BCDict[key1+key].second;
                                        };

                                                
                                };
                        };

                };

                if(globalBCDict[key]=="OUTLET"){

                        wp rho=params["RHO"];
                        wp mDot;
                        wp uFree,vFree;

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){
                                        
                                        uFree=vars["u"](i,j);
                                        vFree=vars["v"](i,j);

                                        if(key=="BCW" || key=="BCE"){
                                              mDot=rho*sig[key]*(mesh->normalsI(i+inc[key].first.first,j+inc[key].first.second).
                                                                         dotProduct(Vec2(uFree,vFree)));
                                        }
                                        else if(key=="BCS" || key=="BCN"){
                                              mDot=rho*sig[key]*(mesh->normalsJ(i+inc[key].first.first,j+inc[key].first.second).
                                                                         dotProduct(Vec2(uFree,vFree)));
                                        };

                                        for(auto key1:vars.getKeys()){
                                                 if(key=="p") continue;
                                                 changeMatElement(i,j,i,j,mDot,coeffs[key1]);
                                        };
                                                
                                };
                        };

                };


        };

};

void icoNSSolver::QDot(){

};

void icoNSSolver::lSolve(){
        SIMPLEDriver();
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

};

void icoNSSolver::getResidual(){

};

bool icoNSSolver::isConverged(){
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
        addPressureSource();

        //Apply BC
        applyBC();

        //Solve x momentum
        ls->setMatFunc(coeffs["u"]);
        ls->setRHSFunc(RHS["u"]);
        ls->solve(vars["u"]);

        //Solve y momentum
        ls->setMatFunc(coeffs["v"]);
        ls->setRHSFunc(RHS["v"]);
        ls->solve(vars["v"]);

        RhieChow();
        assemblePressureCorrectionEqn();

        //Solve pressure-correction eqn
        vars["p'"].initVal(0.0);
        ls->setMatFunc(coeffs["p'"]);
        ls->setRHSFunc(RHS["p'"]);
        ls->solve(vars["p'"]);

        correctPressure();
        correctMassFlux();
        correctVelocity();

//        std::ofstream outfile("out.dat");
//        for(int j=1;j<=vars["p"].intVect[1];j++){
//                for(int i=1;i<=vars["p"].intVect[0];i++){
//                        outfile<<mesh->cc(i,j).x<<"\t"<<mesh->cc(i,j).y<<"\t"<<vars["p"](i,j)<<std::endl; 
//                 };
//                 outfile<<std::endl;
//        };
//
//        outfile.close();
};      

void icoNSSolver::computeDiffusiveFlux(std::string varName){

        Vec2 grad;
        wp d;
        Node pL,pR;
        Vec2 vCC,vCD;
        wp E,T;
        wp gradCC,gradCD;

        wp mu=params["MU"];

        for(int j=1;j<=IFaceMDot.intVect[1];j++){
                for(int i=2;i<IFaceMDot.intVect[0];i++){
                        pL=mesh->cc(i-1,j);
                        pR=mesh->cc(i,j);
                        vCC=Vec2(pR.x-pL.x,pR.y-pL.y);
                        d=vCC.norm2();
                        E=mesh->normalsI(i,j).norm2();
                        vCD=mesh->normalsI(i,j)-vCC*E;
                        T=vCD.norm2();
                        grad=scheme.adjustedFaceGradient(vars[varName],*mesh,i,j,'x');
                        gradCC=(vars[varName](i,j)-vars[varName](i-1,j))/d;
                        gradCD=grad.dotProduct(vCD);

                        //Update cell coefficients of i-1,j
                        changeMatElement(i-1,j,i-1,j,gradCC*mu*E,coeffs[varName]);
                        changeMatElement(i-1,j,i,j,-gradCC*mu*E,coeffs[varName]);
                        RHS[varName](i-1,j)+=gradCD*mu;
                        
                        //Update cell coefficients of i,j
                        changeMatElement(i,j,i,j,gradCC*mu*E,coeffs[varName]);
                        changeMatElement(i,j,i-1,j,-gradCC*mu*E,coeffs[varName]);
                        RHS[varName](i,j)-=gradCD*mu;

                };
        };


        for(int j=2;j<JFaceMDot.intVect[1];j++){
                for(int i=1;i<=JFaceMDot.intVect[0];i++){
                        pL=mesh->cc(i,j-1);
                        pR=mesh->cc(i,j);
                        vCC=Vec2(pR.x-pL.x,pR.y-pL.y);
                        d=vCC.norm2();
                        E=mesh->normalsJ(i,j).norm2();
                        vCD=mesh->normalsJ(i,j)-vCC*E;
                        T=vCD.norm2();
                        grad=scheme.adjustedFaceGradient(vars[varName],*mesh,i,j,'x');
                        gradCC=(vars[varName](i,j)-vars[varName](i,j-1))/d;
                        gradCD=grad.dotProduct(vCD);

                        //Update cell coefficients of i,j-1
                        changeMatElement(i,j-1,i,j-1,gradCC*mu*E,coeffs[varName]);
                        changeMatElement(i,j-1,i,j,-gradCC*mu*E,coeffs[varName]);
                        RHS[varName](i,j-1)+=gradCD*mu;
                        
                        //Update cell coefficients of i,j
                        changeMatElement(i,j,i,j,gradCC*mu*E,coeffs[varName]);
                        changeMatElement(i,j,i,j-1,-gradCC*mu*E,coeffs[varName]);
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
                        for(auto key:coeffs.getKeys()){
                               if(key=="p") continue;
                               for(auto r:recon){ 
                                        changeMatElement(i-1,j,i+r.first,j,IFaceMDot(i,j)*r.second,coeffs[key]);
                                        changeMatElement(i,j,i+r.first,j,-IFaceMDot(i,j)*r.second,coeffs[key]);
                               };
                        };


                        if(scheme.getScheme(faceReconstruction)==QUICK && i>2 && i<IFaceMDot.intVect[0]-1){
                              
                              for(auto key:coeffs.getKeys()){
                                        if(key=="p") continue; 
                                        for(auto r:recon){       
                                                RHS[key](i-1,j)+=vars[key](i+r.first,j)*IFaceMDot(i,j)*r.second;
                                                RHS[key](i,j)-=vars[key](i+r.first,j)*IFaceMDot(i,j)*r.second;
                                        };

                              };

                              recon=scheme.QUICKReconstruction(IFaceMDot(i,j));

                              for(auto key:coeffs.getKeys()){
                                        if(key=="p") continue; 
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
                        for(auto key:coeffs.getKeys()){
                               if(key=="p") continue;
                               for(auto r:recon){ 
                                        changeMatElement(i,j-1,i,j+r.first,JFaceMDot(i,j)*r.second,coeffs[key]);
                                        changeMatElement(i,j,i,j+r.first,-JFaceMDot(i,j)*r.second,coeffs[key]);
                               };
                        };


                        if(scheme.getScheme(faceReconstruction)==QUICK && j>2 && j<JFaceMDot.intVect[1]-1){
                              
                              for(auto key:coeffs.getKeys()){
                                        if(key=="p") continue; 
                                        for(auto r:recon){       
                                                RHS[key](i,j-1)+=vars[key](i,j+r.first)*JFaceMDot(i,j)*r.second;
                                                RHS[key](i,j)-=vars[key](i,j+r.first)*JFaceMDot(i,j)*r.second;
                                        };

                              };

                              recon=scheme.QUICKReconstruction(JFaceMDot(i,j));
                              
                              for(auto key:coeffs.getKeys()){
                                        if(key=="p") continue;          
                                        for(auto r:recon){       
                                                RHS[key](i,j-1)-=vars[key](i,j+r.first)*JFaceMDot(i,j)*r.second;
                                                RHS[key](i,j)+=vars[key](i,j+r.first)*JFaceMDot(i,j)*r.second;
                                        };

                              };

                        };

                };
        };

};

void icoNSSolver::computeGradP(){

        Vec2 gradP;
        for(int j=1;j<=vars["gradP"].intVect[1];j++){
                for(int i=1;i<=vars["gradP"].intVect[0];i++){
                        gradP=scheme.greenGaussCellBased(vars["p"],*mesh,i,j);
                        vars["gradPX"](i,j)=gradP.dotProduct(Vec2(1.0,0.0));
                        vars["gradPY"](i,j)=gradP.dotProduct(Vec2(0.0,1.0));
                };
        };

};

void icoNSSolver::addPressureSource(){

        Vec2 gradP;
        wp gradPX,gradPY;

        for(int j=2;j<vars["u"].intVect[1];j++){
                for(int i=2;i<vars["u"].intVect[0];i++){
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
       
        wp Ef;
        Vec2 dCC;
        Node ccL,ccR;
        wp DLU,DRU,DLV,DRV;

        for(int j=1;j<=IFaceMDot.intVect[1];j++){
                for(int i=2;i<IFaceMDot.intVect[0];i++){
                        ccL=mesh->cc(i-1,j);
                        ccR=mesh->cc(i,j);
                        dCC=Vec2(ccR.x-ccL.x,ccR.y-ccL.y);
                        DLU=mesh->volumes(i-1,j)/findMatElement(i-1,j,i-1,j,coeffs["u"]);
                        DRU=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["u"]);
                        DLV=mesh->volumes(i-1,j)/findMatElement(i-1,j,i-1,j,coeffs["v"]);
                        DRV=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["v"]);

                        Ef=dCC.en.x*scheme.averageValue<wp>(DLU,DRU,*mesh,i,j,'x')*mesh->normalsI(i,j).en.x+
                           dCC.en.y*scheme.averageValue<wp>(DLV,DRV,*mesh,i,j,'x')*mesh->normalsI(i,j).en.y;

                        Ef=Ef/pow(dCC.norm2(),2.0);
                        changeMatElement(i-1,j,i-1,j,-Ef,coeffs["p'"]);
                        changeMatElement(i-1,j,i,j,Ef,coeffs["p'"]);

                        changeMatElement(i,j,i,j,-Ef,coeffs["p'"]);
                        changeMatElement(i,j,i-1,j,Ef,coeffs["p'"]);

                };
        };

        for(int j=2;j<JFaceMDot.intVect[1];j++){
                for(int i=1;i<=JFaceMDot.intVect[0];i++){
                        ccL=mesh->cc(i,j-1);
                        ccR=mesh->cc(i,j);
                        dCC=Vec2(ccR.x-ccL.x,ccR.y-ccL.y);
                        DLU=mesh->volumes(i,j-1)/findMatElement(i,j-1,i,j-1,coeffs["u"]);
                        DRU=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["u"]);
                        DLV=mesh->volumes(i,j-1)/findMatElement(i,j-1,i,j-1,coeffs["v"]);
                        DRV=mesh->volumes(i,j)/findMatElement(i,j,i,j,coeffs["v"]);

                        Ef=dCC.en.x*scheme.averageValue<wp>(DLU,DRU,*mesh,i,j,'y')*mesh->normalsJ(i,j).en.x+
                           dCC.en.y*scheme.averageValue<wp>(DLV,DRV,*mesh,i,j,'y')*mesh->normalsJ(i,j).en.y;

                        Ef=Ef/pow(dCC.norm2(),2.0);
                        changeMatElement(i,j-1,i,j-1,-Ef,coeffs["p'"]);
                        changeMatElement(i,j-1,i,j,Ef,coeffs["p'"]);

                        changeMatElement(i,j,i,j,-Ef,coeffs["p'"]);
                        changeMatElement(i,j,i,j-1,Ef,coeffs["p'"]);

                };
        };

        for(int j=2;j<vars["p'"].intVect[1];j++){
                for(int i=2;i<vars["p'"].intVect[0];i++){

                        RHS["p'"](i,j)+=IFaceMDot(i+1,j)-IFaceMDot(i,j)
                                        +JFaceMDot(i,j+1)-JFaceMDot(i,j);
                };
        };

};

void icoNSSolver::correctPressure(){

};

void icoNSSolver::correctMassFlux(){

};

void icoNSSolver::correctVelocity(){


};

