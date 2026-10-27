#include"comNSSolver.h"

comNSSolver::comNSSolver(){
        std::ofstream outfile("res.dat");
        outfile<<"Continuity\tXMom\tYMom\tEnergy"<<std::endl;
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
        (*this).setVar("rho",{CELL,CELL});
        (*this).setVar("u",{CELL,CELL});
        (*this).setVar("v",{CELL,CELL});
        (*this).setVar("p",{CELL,CELL});

        for(auto key:vars.getKeys()){
                x_qp_left[key]=Grid<wp>(mesh->box,{NODE,CELL});
                x_qp_right[key]=Grid<wp>(mesh->box,{NODE,CELL});
                y_qp_left[key]=Grid<wp>(mesh->box,{CELL,NODE});
                y_qp_right[key]=Grid<wp>(mesh->box,{CELL,NODE});
                F[key]=Grid<wp>(mesh->box,{NODE,CELL});
                G[key]=Grid<wp>(mesh->box,{CELL,NODE});
        };

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
                                };
                                continue;
               }
               else if(data=="BOUNDARYCONDITIONS"){
                                logger.log("Reading BC",1);
                                while(infile1>>data){
                                        if(data.front()=='#') break;
                                        infile1>>data1;
                                        if(data1=="SLIPWALL"){
                                                globalBCDict[data]=data1;
                                                infile1>>dataNum;
                                                scheme.BCDict["rho"+data]={neumann,0.0};
                                                infile1>>dataNum;
                                                scheme.BCDict["u"+data]={dirichlet,dataNum};
                                                infile1>>dataNum;
                                                scheme.BCDict["v"+data]={dirichlet,dataNum};
                                                infile1>>dataNum;
                                                scheme.BCDict["p"+data]={neumann,0.0};
                                        }
                                        if(data1=="NOSLIPWALL"){
                                                globalBCDict[data]=data1;
                                                infile1>>dataNum;
                                                scheme.BCDict["rho"+data]={neumann,0.0};
                                                infile1>>dataNum;
                                                scheme.BCDict["u"+data]={dirichlet,dataNum};
                                                infile1>>dataNum;
                                                scheme.BCDict["v"+data]={dirichlet,dataNum};
                                                infile1>>dataNum;
                                                scheme.BCDict["p"+data]={neumann,0.0};
                                        }
                                        else if(data1=="SUPERSONICINLET"){
                                                globalBCDict[data]=data1;
                                                infile1>>dataNum;
                                                scheme.BCDict["rho"+data]={dirichlet,dataNum};
                                                infile1>>dataNum;
                                                scheme.BCDict["u"+data]={dirichlet,dataNum};
                                                infile1>>dataNum;
                                                scheme.BCDict["v"+data]={dirichlet,dataNum};
                                                infile1>>dataNum;
                                                scheme.BCDict["p"+data]={dirichlet,dataNum};
                                        }
                                        else if(data1=="SUPERSONICOUTLET"){
                                                globalBCDict[data]=data1;
                                                infile1>>dataNum;
                                                scheme.BCDict["rho"+data]={extrapolate,0.0};
                                                infile1>>dataNum;
                                                scheme.BCDict["u"+data]={extrapolate,0.0};
                                                infile1>>dataNum;
                                                scheme.BCDict["v"+data]={extrapolate,0.0};
                                                infile1>>dataNum;
                                                scheme.BCDict["p"+data]={extrapolate,0.0};
                                        }
                                        else if(data1=="SUBSONICINLET"){
                                                globalBCDict[data]=data1;
                                                infile1>>dataNum;
                                                scheme.BCDict["rho"+data]={dirichlet,dataNum};
                                                infile1>>dataNum;
                                                scheme.BCDict["u"+data]={dirichlet,dataNum};
                                                infile1>>dataNum;
                                                scheme.BCDict["v"+data]={dirichlet,dataNum};
                                                infile1>>dataNum;
                                                scheme.BCDict["p"+data]={extrapolate,0.0};
                                        }
                                        else if(data1=="SUBSONICOUTLET"){
                                                globalBCDict[data]=data1;
                                                infile1>>dataNum;
                                                scheme.BCDict["rho"+data]={extrapolate,0.0};
                                                infile1>>dataNum;
                                                scheme.BCDict["u"+data]={extrapolate,0.0};
                                                infile1>>dataNum;
                                                scheme.BCDict["v"+data]={extrapolate,0.0};
                                                infile1>>dataNum;
                                                scheme.BCDict["p"+data]={dirichlet,dataNum};
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

        res["rho"]=0.0;
        res["u"]=0.0;
        res["v"]=0.0;
        res["p"]=0.0;

        resScale["rho"]=params["RHOINF"]*Vec2(params["UINF"],params["VINF"]).norm2();
        resScale["u"]=params["RHOINF"]*pow(Vec2(params["UINF"],params["VINF"]).norm2(),2.0);
        resScale["v"]=params["RHOINF"]*pow(Vec2(params["UINF"],params["VINF"]).norm2(),2.0);
        resScale["p"]=0.5*params["RHOINF"]*pow(Vec2(params["UINF"],params["VINF"]).norm2(),3.0)
                      +params["GM"]*params["PINF"]/(params["GM"]-1.0);


};

void comNSSolver::initialCondition(){

        std::ifstream infile(initFile);
        std::string temp;

        infile>>temp;
        for(int j=1;j<=vars["rho"].intVect[1];j++){
                for(int i=1;i<=vars["rho"].intVect[0];i++){
                        infile>>vars["rho"](i,j);
                };
        };

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

void comNSSolver::applyBC(){

        dictionary<std::string,std::pair<int,int>>lo;
        dictionary<std::string,std::pair<int,int>>hi;
        dictionary<std::string,std::pair<int,int>>inc;
        dictionary<std::string,wp>sig;
     
        lo["BCW"]={1,1};
        hi["BCW"]={1,mesh->box.jmx-1};
        inc["BCW"]={-1,0};
        sig["BCW"]=-1.0;

        lo["BCE"]={mesh->box.imx-1,1};
        hi["BCE"]={mesh->box.imx-1,mesh->box.jmx-1};
        inc["BCE"]={1,0};
        sig["BCE"]=1.0;

        lo["BCS"]={1,1};
        hi["BCS"]={mesh->box.imx-1,1};
        inc["BCS"]={0,-1};
        sig["BCS"]=-1.0;

        lo["BCN"]={1,mesh->box.jmx-1};
        hi["BCN"]={mesh->box.imx-1,mesh->box.jmx-1};
        inc["BCN"]={0,1};
        sig["BCN"]=1.0;

        zeroFluxBC["BCS"]=false;
        zeroFluxBC["BCN"]=false;
        zeroFluxBC["BCE"]=false;
        zeroFluxBC["BCW"]=false;

        wp rhoTemp,uTemp,vTemp,pTemp;
        wp velParallel,velNormal;
        wp S;
        Vec2 nHat,pHat;
        Vec2 velNet;

        Grid<wp>&rho=vars["rho"];
        Grid<wp>&u=vars["u"];
        Grid<wp>&v=vars["v"];
        Grid<wp>&p=vars["p"];

        for(auto key:globalBCDict.getKeys()){
                if(globalBCDict[key]=="SLIPWALL"){
                       
                        zeroFluxBC[key]=true;

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){

                                        uTemp=u(i,j);
                                        vTemp=v(i,j);
                                        
                                        if(key=="BCW"){
                                                S=mesh->normalsI(i,j).norm2();
                                                nHat=mesh->normalsI(i,j)/S;
                                                pHat=Vec2((*mesh)(i,j+1).x-(*mesh)(i,j).x,
                                                          (*mesh)(i,j+1).y-(*mesh)(i,j).y);
                                                pHat=pHat/S;
                                        }
                                        else if(key=="BCE"){
                                                S=mesh->normalsI(i+1,j).norm2();
                                                nHat=mesh->normalsI(i+1,j)/S;
                                                pHat=Vec2((*mesh)(i+1,j+1).x-(*mesh)(i+1,j).x,
                                                          (*mesh)(i+1,j+1).y-(*mesh)(i+1,j).y);
                                                pHat=pHat/S;
                                        }
                                        else if(key=="BCS"){
                                                S=mesh->normalsJ(i,j).norm2();
                                                nHat=mesh->normalsJ(i,j)/S;
                                                pHat=Vec2((*mesh)(i+1,j).x-(*mesh)(i,j).x,
                                                          (*mesh)(i+1,j).y-(*mesh)(i,j).y);
                                                pHat=pHat/S;
                                        }
                                        else if(key=="BCN"){
                                                S=mesh->normalsJ(i,j+1).norm2();
                                                nHat=mesh->normalsJ(i,j+1)/S;
                                                pHat=Vec2((*mesh)(i+1,j+1).x-(*mesh)(i,j+1).x,
                                                          (*mesh)(i+1,j+1).y-(*mesh)(i,j+1).y);
                                                pHat=pHat/S;
                                        };

                                        velNormal=Vec2(uTemp,vTemp).dotProduct(nHat);
                                        velParallel=Vec2(uTemp,vTemp).dotProduct(pHat);

                                        velNet=(pHat*velParallel)-(nHat*velNormal);

                                        rho(i+inc[key].first,j+inc[key].second)=rho(i,j);
                                        u(i+inc[key].first,j+inc[key].second)=velNet.en.x;
                                        v(i+inc[key].first,j+inc[key].second)=velNet.en.y;
                                        p(i+inc[key].first,j+inc[key].second)=p(i,j);

                                };
                        };

                };

                if(globalBCDict[key]=="NOSLIPWALL"){

                        zeroFluxBC[key]=true;

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){

                                        rho(i+inc[key].first,j+inc[key].second)=rho(i,j);
                                        u(i+inc[key].first,j+inc[key].second)=-u(i,j);
                                        v(i+inc[key].first,j+inc[key].second)=-v(i,j);
                                        p(i+inc[key].first,j+inc[key].second)=p(i,j);

                                };
                        };

                };

                if(globalBCDict[key]=="SUPERSONICINLET"){
                        
                        wp rhoIn=scheme.BCDict["rho"+key].second;
                        wp uIn=scheme.BCDict["u"+key].second;
                        wp vIn=scheme.BCDict["v"+key].second;
                        wp pIn=scheme.BCDict["p"+key].second;

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){

                                        rho(i+inc[key].first,j+inc[key].second)=rhoIn;
                                        u(i+inc[key].first,j+inc[key].second)=uIn;
                                        v(i+inc[key].first,j+inc[key].second)=vIn;
                                        p(i+inc[key].first,j+inc[key].second)=pIn;

                                };
                        };

                };

                if(globalBCDict[key]=="SUPERSONICOUTLET"){

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){

                                        rho(i+inc[key].first,j+inc[key].second)=rho(i,j);
                                        u(i+inc[key].first,j+inc[key].second)=u(i,j);
                                        v(i+inc[key].first,j+inc[key].second)=v(i,j);
                                        p(i+inc[key].first,j+inc[key].second)=p(i,j);

                                };
                        };

                };

                if(globalBCDict[key]=="SUBSONICINLET"){

                        wp rhoIn=scheme.BCDict["rho"+key].second;
                        wp uIn=scheme.BCDict["u"+key].second;
                        wp vIn=scheme.BCDict["v"+key].second;

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){

                                        rho(i+inc[key].first,j+inc[key].second)=rhoIn;
                                        u(i+inc[key].first,j+inc[key].second)=uIn;
                                        v(i+inc[key].first,j+inc[key].second)=vIn;
                                        p(i+inc[key].first,j+inc[key].second)=p(i,j);

                                };
                        };

                };

                if(globalBCDict[key]=="SUBSONICOUTLET"){

                        wp pIn=scheme.BCDict["p"+key].second;

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){

                                        rho(i+inc[key].first,j+inc[key].second)=rho(i,j);
                                        u(i+inc[key].first,j+inc[key].second)=u(i,j);
                                        v(i+inc[key].first,j+inc[key].second)=v(i,j);
                                        p(i+inc[key].first,j+inc[key].second)=pIn;

                                };
                        };

                };

        };

};

void comNSSolver::QDot(){
       
       applyBC();
       reconstructStates();
       computeFluxes();
       getResidual();
        
};

void comNSSolver::lSolve(){

};

void comNSSolver::computeTimeStep(Grid<wp>&dt){

        Vec2 nHat;
        wp S;
        wp contraVel,c;
        wp gm=params["GM"];

        Grid<wp>&rho=vars["rho"];
        Grid<wp>&u=vars["u"];
        Grid<wp>&v=vars["v"];
        Grid<wp>&p=vars["p"];

        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        
                        dt(i,j)=0.0;
                        
                        //Left face
                        S=mesh->normalsI(i,j).norm2();
                        nHat=mesh->normalsI(i,j)/S;
                        contraVel=fabs(nHat.dotProduct((Vec2(u(i-1,j),v(i-1,j))+
                                                   Vec2(u(i,j),v(i,j)))*0.5));
                        c=0.5*(sqrt(gm*p(i-1,j)/rho(i-1,j))
                              +sqrt(gm*p(i,j)/rho(i,j)));

                        dt(i,j)=dt(i,j)+(contraVel+c)*S;

                        //Right face
                        S=mesh->normalsI(i+1,j).norm2();
                        nHat=mesh->normalsI(i+1,j)/S;
                        contraVel=fabs(nHat.dotProduct((Vec2(u(i,j),v(i,j))+
                                                   Vec2(u(i+1,j),v(i+1,j)))*0.5));
                        c=0.5*(sqrt(gm*p(i,j)/rho(i,j))
                              +sqrt(gm*p(i+1,j)/rho(i+1,j)));

                        dt(i,j)=dt(i,j)+(contraVel+c)*S;

                        //Bottom face
                        S=mesh->normalsJ(i,j).norm2();
                        nHat=mesh->normalsJ(i,j)/S;
                        contraVel=fabs(nHat.dotProduct((Vec2(u(i,j-1),v(i,j-1))+
                                                   Vec2(u(i,j),v(i,j)))*0.5));
                        c=0.5*(sqrt(gm*p(i,j-1)/rho(i,j-1))
                              +sqrt(gm*p(i,j)/rho(i,j)));

                        dt(i,j)=dt(i,j)+(contraVel+c)*S;

                        //Top face
                        S=mesh->normalsJ(i,j+1).norm2();
                        nHat=mesh->normalsJ(i,j+1)/S;
                        contraVel=fabs(nHat.dotProduct((Vec2(u(i,j),v(i,j))+
                                                   Vec2(u(i,j+1),v(i,j+1)))*0.5));
                        c=0.5*(sqrt(gm*p(i,j)/rho(i,j))
                              +sqrt(gm*p(i,j+1)/rho(i,j+1)));

                        dt(i,j)=dt(i,j)+(contraVel+c)*S;

                        dt(i,j)=mesh->volumes(i,j)/dt(i,j);

                };
        };

        if(scheme.getScheme(timeSteppingLocality)==GLOBAL){
                wp mindt=BIG;

                for(int j=1;j<=vars["u"].intVect[1];j++){
                        for(int i=1;i<=vars["u"].intVect[0];i++){
                             
                             if(dt(i,j)<mindt) mindt=dt(i,j);   

                        };
                };

                for(int j=1;j<=vars["u"].intVect[1];j++){
                        for(int i=1;i<=vars["u"].intVect[0];i++){
                             
                                dt(i,j)=mindt;

                        };
                };

        };

};

void comNSSolver::updateVars(Grid<wp>&dt, wp storeFactor){
        
        std::vector<wp>qp;
        std::vector<wp>qc;

        std::vector<std::string>keys=vars.getKeys();

        Grid<wp>&rho=vars["rho"];
        Grid<wp>&u=vars["u"];
        Grid<wp>&v=vars["v"];
        Grid<wp>&p=vars["p"];

        wp gm=params["GM"];
        wp RGas=params["RGAS"];

        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        qp={rho(i,j),u(i,j),v(i,j),p(i,j)};
                        qc=primitiveToConservative(qp,gm,RGas);

                        rho(i,j)=qc[0];
                        u(i,j)=qc[1];
                        v(i,j)=qc[2];
                        p(i,j)=qc[3];

                };
        };

        for(auto key:keys){
                for(int j=1;j<=vars["u"].intVect[1];j++){
                        for(int i=1;i<=vars["u"].intVect[0];i++){
                                vars[key](i,j)=vars[key](i,j)+varsDot[key](i,j)*(dt(i,j)*storeFactor);
                        };

                };
        };

        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        qc={rho(i,j),u(i,j),v(i,j),p(i,j)};
                        qp=conservativeToPrimitive(qc,gm,RGas);

                        rho(i,j)=qp[0];
                        u(i,j)=qp[1];
                        v(i,j)=qp[2];
                        p(i,j)=qp[3];

                };
        };

};

void comNSSolver::updateVars(Grid<wp>&dt, wp storeFactor, std::vector<Grid<wp>>& uStore){
        std::vector<wp>qp;
        std::vector<wp>qc;

        std::vector<std::string>keys=vars.getKeys();

        Grid<wp>&rho=vars["rho"];
        Grid<wp>&u=vars["u"];
        Grid<wp>&v=vars["v"];
        Grid<wp>&p=vars["p"];

        wp gm=params["GM"];
        wp RGas=params["RGAS"];

        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        qp={uStore[0](i,j),uStore[1](i,j),uStore[2](i,j),uStore[3](i,j)};
                        qc=primitiveToConservative(qp,gm,RGas);

                        rho(i,j)=qc[0];
                        u(i,j)=qc[1];
                        v(i,j)=qc[2];
                        p(i,j)=qc[3];

                };
        };

        for(auto key:keys){
                for(int j=1;j<=vars["u"].intVect[1];j++){
                        for(int i=1;i<=vars["u"].intVect[0];i++){
                                vars[key](i,j)=vars[key](i,j)+varsDot[key](i,j)*(dt(i,j)*storeFactor);
                        };

                };
        };

        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        qc={rho(i,j),u(i,j),v(i,j),p(i,j)};
                        qp=conservativeToPrimitive(qc,gm,RGas);

                        rho(i,j)=qp[0];
                        u(i,j)=qp[1];
                        v(i,j)=qp[2];
                        p(i,j)=qp[3];

                };
        };

};

void comNSSolver::updateVars(Grid<wp>&dt, wp storeFactor, std::vector<Grid<wp>>& uStore, std::vector<Grid<wp>>& rStore){
        std::vector<wp>qp;
        std::vector<wp>qc;

        std::vector<std::string>keys=vars.getKeys();

        Grid<wp>&rho=vars["rho"];
        Grid<wp>&u=vars["u"];
        Grid<wp>&v=vars["v"];
        Grid<wp>&p=vars["p"];

        wp gm=params["GM"];
        wp RGas=params["RGAS"];

        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        qp={uStore[0](i,j),uStore[1](i,j),uStore[2](i,j),uStore[3](i,j)};
                        qc=primitiveToConservative(qp,gm,RGas);

                        rho(i,j)=qc[0];
                        u(i,j)=qc[1];
                        v(i,j)=qc[2];
                        p(i,j)=qc[3];

                };
        };

        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){

                        for(int n=0;n<keys.size();n++){
                                vars[keys[n]](i,j)=vars[keys[n]](i,j)+rStore[n](i,j)*(dt(i,j)*storeFactor);
                        };

                };
        };

        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        qc={rho(i,j),u(i,j),v(i,j),p(i,j)};
                        qp=conservativeToPrimitive(qc,gm,RGas);

                        rho(i,j)=qp[0];
                        u(i,j)=qp[1];
                        v(i,j)=qp[2];
                        p(i,j)=qp[3];

                };
        };

};

wp comNSSolver::getResNorm(){
        
        std::vector<std::string>keys=vars.getKeys();

        for(auto key:keys){
                res[key]=0.0;
        };

        for(auto key:keys){
                wp&resLocal=res[key];
                Grid<wp>&varsDotLocal=varsDot[key];
                for(int j=1;j<=vars["u"].intVect[1];j++){
                        for(int i=1;i<=vars["u"].intVect[0];i++){
                                resLocal=resLocal+pow(varsDotLocal(i,j)*mesh->volumes(i,j),2.0);
                        };
                };
        };

        for(auto key:keys){
                res[key]=sqrt(res[key])/resScale[key];
          };

        std::ofstream outfile("res.dat",std::ios::app);

        for(auto key:keys){
                outfile<<res[key]<<"\t";
        };

        outfile<<std::endl;
        outfile.close();

};

void comNSSolver::getResidual(){

        std::vector<std::string>keys=vars.getKeys();

        for(auto key:keys){
                Grid<wp>&varsDotLocal=varsDot[key];
                Grid<wp>&FLocal=F[key];
                Grid<wp>&GLocal=G[key];
                for(int j=1;j<=vars["u"].intVect[1];j++){
                        for(int i=1;i<=vars["u"].intVect[0];i++){
                                varsDotLocal(i,j)=-(FLocal(i+1,j)-FLocal(i,j)+GLocal(i,j+1)-GLocal(i,j))/mesh->volumes(i,j);
                        };
                };
        };

};

bool comNSSolver::isConverged(){
        std::vector<std::string>keys=vars.getKeys();

        wp conv=0.0;

        for(auto key:keys){
                conv=conv+pow(res[key],2.0);
        };

        return sqrt(conv)<params["TOL"];

};

void comNSSolver::reconstructStates(){
        
        std::vector<std::string>keys=vars.getKeys();

        if(scheme.getScheme(stateReconstruction)==NONE){
                
                for(int j=1;j<=vars["u"].intVect[1];j++){
                        for(int i=1;i<=vars["u"].intVect[0]+1;i++){
                                for(auto key:keys){ 
                                        scheme.defaultReconstruction(vars[key],x_qp_left[key],x_qp_right[key],i,j,{1,0});
                                };
                        };
                };

                for(int j=1;j<=vars["u"].intVect[1]+1;j++){
                        for(int i=1;i<=vars["u"].intVect[0];i++){
                                for(auto key:keys){ 
                                        scheme.defaultReconstruction(vars[key],y_qp_left[key],y_qp_right[key],i,j,{0,1});
                                };
                        };
                };

        };

};

void comNSSolver::computeFluxes(){

        if(scheme.getScheme(fluxCalc)==AUSM){
                
                scheme.ausmFlux(params,{zeroFluxBC["BCW"],zeroFluxBC["BCE"]},mesh->normalsI,F,x_qp_left,x_qp_right,{1,0});
                scheme.ausmFlux(params,{zeroFluxBC["BCS"],zeroFluxBC["BCN"]},mesh->normalsJ,G,y_qp_left,y_qp_right,{0,1});
        
        };

};

std::vector<wp> comNSSolver::conservativeToPrimitive(std::vector<wp>qc, wp gm, wp RGas){
        std::vector<wp>qp;
        qp.resize(4);

        wp cvT;

        qp[0]=qc[0];
        qp[1]=qc[1]/qc[0];
        qp[2]=qc[2]/qc[0];
        cvT=(qc[3]/qc[0])-0.5*(qp[1]*qp[1]+qp[2]*qp[2]);
        qp[3]=qp[0]*cvT*(gm-1.0);

        return qp;

};

std::vector<wp> comNSSolver::primitiveToConservative(std::vector<wp>qp, wp gm, wp RGas){
        std::vector<wp>qc;
        qc.resize(4);

        wp T;

        qc[0]=qp[0];
        qc[1]=qp[0]*qp[1];
        qc[2]=qp[0]*qp[2];
        T=qp[3]/(RGas*qp[0]);
        qc[3]=qp[0]*((RGas/(gm-1.0))*T+0.5*(qp[1]*qp[1]+qp[2]*qp[2]));

        return qc;

};

