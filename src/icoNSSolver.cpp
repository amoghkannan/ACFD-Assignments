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
        RHS["u"]=vars["u"];
        RHS["v"]=vars["v"];
        RHS["p"]=vars["p"];
        coeffs["u"]=Grid<boundMatRow>(Box(imx,jmx,1,1,1,1),{CELL,CELL});
        coeffs["v"]=Grid<boundMatRow>(Box(imx,jmx,1,1,1,1),{CELL,CELL});
        coeffs["p"]=Grid<boundMatRow>(Box(imx,jmx,1,1,1,1),{CELL,CELL});
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
};

void icoNSSolver::applyBC(){

        dictionary<std::string,std::pair<int,int>>lo;
        dictionary<std::string,std::pair<int,int>>hi;
        dictionary<std::string,std::pair<std::pair<int,int>,std::pair<int,int>>>inc;
        dictionary<std::string,wp>sign;
     
        lo["BCW"]={1,1};
        hi["BCW"]={1,mesh->box.jmx-1};
        inc["BCW"]={{0,0},{0,1}};
        sign["BCW"]=-1.0;

        lo["BCE"]={mesh->box.imx-1,1};
        hi["BCE"]={mesh->box.imx-1,mesh->box.jmx-1};
        inc["BCE"]={{1,0},{1,1}};
        sign["BCE"]=1.0;

        lo["BCS"]={1,1};
        hi["BCS"]={mesh->box.imx-1,1};
        inc["BCS"]={{0,0},{1,0}};
        sign["BCS"]=-1.0;

        lo["BCN"]={1,mesh->box.jmx-1};
        hi["BCN"]={mesh->box.imx-1,mesh->box.jmx-1};
        inc["BCN"]={{0,1},{1,1}};
        sign["BCN"]=1.0;

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

                };
        };


        //Compute coefficients
        computeDiffusiveFlux("u");
        computeDiffusiveFlux("v");

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

};

void icoNSSolver::computeDiffusiveFlux(std::string varName){

        wp gradL,gradR,grad;
        wp d;
        Vec2 dCC;
        wp SCC,SCD;

        wp mu=params["MU"];

        for(int j=1;j<=IFaceMDot.intVect[1];j++){
                for(int i=2;i<IFaceMDot.intVect[0];i++){
                        dCC=Vec2(mesh->cc(i-1,j),mesh->cc(i,j));
                        d=dCC.norm2();
                        SCC=mesh->normalsI(i,j).dotProduct(dCC)/d;
                        SCD=(mesh->normalsI(i,j)-(dCC*SCC)/d).norm2();
                        gradL=scheme.greenGaussCellBased(vars[varName],*mesh,dCC,i-1,j);
                        gradR=scheme.greenGaussCellBased(vars[varName],*mesh,dCC,i,j);
        
                        //Update cell coefficients of i-1,j
                        changeMatElement(i-1,j,i-1,j,SCC*mu/d,coeffs[varName]);
                        changeMatElement(i-1,j,i,j,-SCC*mu/d,coeffs[varName]);
                        RHS[varName](i-1,j)+=0.5*(gradL+gradR)*mu*SCD;
                        
                        //Update cell coefficients of i,j
                        changeMatElement(i,j,i,j,SCC*mu/d,coeffs[varName]);
                        changeMatElement(i,j,i-1,j,-SCC*mu/d,coeffs[varName]);
                        RHS[varName](i,j)-=0.5*(gradL+gradR)*mu*SCD;

                };
        };


        for(int j=2;j<JFaceMDot.intVect[1];j++){
                for(int i=1;i<=JFaceMDot.intVect[0];i++){
                        dCC=Vec2(mesh->cc(i,j-1),mesh->cc(i,j));
                        d=dCC.norm2();
                        SCC=mesh->normalsJ(i,j).dotProduct(dCC)/d;
                        SCD=(mesh->normalsJ(i,j)-(dCC*SCC)/d).norm2();
                        gradL=scheme.greenGaussCellBased(vars[varName],*mesh,dCC,i,j-1);
                        gradR=scheme.greenGaussCellBased(vars[varName],*mesh,dCC,i,j);
        
                        //Update cell coefficients of i,j-1
                        changeMatElement(i,j-1,i,j-1,SCC*mu/d,coeffs[varName]);
                        changeMatElement(i,j-1,i,j,-SCC*mu/d,coeffs[varName]);
                        RHS[varName](i,j-1)+=0.5*(gradL+gradR)*mu*SCD;
                        
                        //Update cell coefficients of i,j
                        changeMatElement(i,j,i,j,SCC*mu/d,coeffs[varName]);
                        changeMatElement(i,j,i,j-1,-SCC*mu/d,coeffs[varName]);
                        RHS[varName](i,j)-=0.5*(gradL+gradR)*mu*SCD;

                };
        };

};

