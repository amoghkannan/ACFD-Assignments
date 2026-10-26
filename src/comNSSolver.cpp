#include"comNSSolver.h"

comNSSolver::comNSSolver(){
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
        (*this).setVar("rho",{CELL,CELL});
        (*this).setVar("u",{CELL,CELL});
        (*this).setVar("v",{CELL,CELL});
        (*this).setVar("p",{CELL,CELL});

        for(auto key:vars.getKeys()){
                x_qp_left[key]=Grid<wp>(mesh->box,{NODE,CELL});
                x_qp_right[key]=Grid<wp>(mesh->box,{NODE,CELL});
                y_qp_left[key]=Grid<wp>(mesh->box,{CELL,NODE});
                y_qp_right[key]=Grid<wp>(mesh->box,{CELL,NODE});
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

        wp rhoTemp,uTemp,vTemp,pTemp;
        wp velParallel,velNormal;
        wp S;
        Vec2 nHat,pHat;
        Vec2 velNet;

        for(auto key:globalBCDict.getKeys()){
                if(globalBCDict[key]=="SLIPWALL"){
                        
                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){

                                        uTemp=vars["u"](i,j);
                                        vTemp=vars["v"](i,j);
                                        
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

                                        vars["rho"](i+inc[key].first,j+inc[key].second)=vars["rho"](i,j);
                                        vars["u"](i+inc[key].first,j+inc[key].second)=velNet.en.x;
                                        vars["v"](i+inc[key].first,j+inc[key].second)=velNet.en.y;
                                        vars["p"](i+inc[key].first,j+inc[key].second)=vars["p"](i,j);

                                };
                        };

                };

                if(globalBCDict[key]=="NOSLIPWALL"){

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){

                                        vars["rho"](i+inc[key].first,j+inc[key].second)=vars["rho"](i,j);
                                        vars["u"](i+inc[key].first,j+inc[key].second)=-vars["u"](i,j);
                                        vars["v"](i+inc[key].first,j+inc[key].second)=-vars["v"](i,j);
                                        vars["p"](i+inc[key].first,j+inc[key].second)=vars["p"](i,j);

                                };
                        };

                };

                if(globalBCDict[key]=="SUPERSONICINLET"){

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){

                                        vars["rho"](i+inc[key].first,j+inc[key].second)=scheme.BCDict["rho"+key].second;
                                        vars["u"](i+inc[key].first,j+inc[key].second)=scheme.BCDict["u"+key].second;
                                        vars["v"](i+inc[key].first,j+inc[key].second)=scheme.BCDict["v"+key].second;
                                        vars["p"](i+inc[key].first,j+inc[key].second)=scheme.BCDict["p"+key].second;

                                };
                        };

                };

                if(globalBCDict[key]=="SUPERSONICOUTLET"){

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){

                                        vars["rho"](i+inc[key].first,j+inc[key].second)=vars["rho"](i,j);
                                        vars["u"](i+inc[key].first,j+inc[key].second)=vars["u"](i,j);
                                        vars["v"](i+inc[key].first,j+inc[key].second)=vars["v"](i,j);
                                        vars["p"](i+inc[key].first,j+inc[key].second)=vars["p"](i,j);

                                };
                        };

                };

                if(globalBCDict[key]=="SUBSONICINLET"){

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){

                                        vars["rho"](i+inc[key].first,j+inc[key].second)=scheme.BCDict["rho"+key].second;
                                        vars["u"](i+inc[key].first,j+inc[key].second)=scheme.BCDict["u"+key].second;
                                        vars["v"](i+inc[key].first,j+inc[key].second)=scheme.BCDict["v"+key].second;
                                        vars["p"](i+inc[key].first,j+inc[key].second)=vars["p"](i,j);

                                };
                        };

                };

                if(globalBCDict[key]=="SUBSONICOUTLET"){

                        for(int j=lo[key].second;j<=hi[key].second;j++){
                                for(int i=lo[key].first;i<=hi[key].first;i++){

                                        vars["rho"](i+inc[key].first,j+inc[key].second)=vars["rho"](i,j);
                                        vars["u"](i+inc[key].first,j+inc[key].second)=vars["u"](i,j);
                                        vars["v"](i+inc[key].first,j+inc[key].second)=vars["v"](i,j);
                                        vars["p"](i+inc[key].first,j+inc[key].second)=scheme.BCDict["p"+key].second;

                                };
                        };

                };

        };

};

void comNSSolver::QDot(){

};

void comNSSolver::lSolve(){

};

void comNSSolver::computeTimeStep(Grid<wp>&dt){

        Vec2 nHat;
        wp S;
        wp contraVel,c;

        for(int j=1;j<=vars["u"].intVect[1];j++){
                for(int i=1;i<=vars["u"].intVect[0];i++){
                        
                        dt(i,j)=0.0;
                        
                        //Left face
                        S=mesh->normalsI(i,j).norm2();
                        nHat=mesh->normalsI(i,j)/S;
                        contraVel=fabs(nHat.dotProduct((Vec2(vars["u"](i-1,j),vars["v"](i-1,j))+
                                                   Vec2(vars["u"](i,j),vars["v"](i,j)))*0.5));
                        c=0.5*(sqrt(params["GM"]*vars["p"](i-1,j)/vars["rho"](i-1,j))
                              +sqrt(params["GM"]*vars["p"](i,j)/vars["rho"](i,j)));

                        dt(i,j)=dt(i,j)+(contraVel+c)*S;

                        //Right face
                        S=mesh->normalsI(i+1,j).norm2();
                        nHat=mesh->normalsI(i+1,j)/S;
                        contraVel=fabs(nHat.dotProduct((Vec2(vars["u"](i,j),vars["v"](i,j))+
                                                   Vec2(vars["u"](i+1,j),vars["v"](i+1,j)))*0.5));
                        c=0.5*(sqrt(params["GM"]*vars["p"](i,j)/vars["rho"](i,j))
                              +sqrt(params["GM"]*vars["p"](i+1,j)/vars["rho"](i+1,j)));

                        dt(i,j)=dt(i,j)+(contraVel+c)*S;

                        //Bottom face
                        S=mesh->normalsJ(i,j).norm2();
                        nHat=mesh->normalsJ(i,j)/S;
                        contraVel=fabs(nHat.dotProduct((Vec2(vars["u"](i,j-1),vars["v"](i,j-1))+
                                                   Vec2(vars["u"](i,j),vars["v"](i,j)))*0.5));
                        c=0.5*(sqrt(params["GM"]*vars["p"](i,j-1)/vars["rho"](i,j-1))
                              +sqrt(params["GM"]*vars["p"](i,j)/vars["rho"](i,j)));

                        dt(i,j)=dt(i,j)+(contraVel+c)*S;

                        //Top face
                        S=mesh->normalsJ(i,j+1).norm2();
                        nHat=mesh->normalsJ(i,j+1)/S;
                        contraVel=fabs(nHat.dotProduct((Vec2(vars["u"](i,j),vars["v"](i,j))+
                                                   Vec2(vars["u"](i,j+1),vars["v"](i,j+1)))*0.5));
                        c=0.5*(sqrt(params["GM"]*vars["p"](i,j)/vars["rho"](i,j))
                              +sqrt(params["GM"]*vars["p"](i,j+1)/vars["rho"](i,j+1)));

                        dt(i,j)=dt(i,j)+(contraVel+c)*S;

                        dt(i,j)=params["CFL"]*mesh->volumes(i,j)/dt(i,j);

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

};

void comNSSolver::updateVars(Grid<wp>&dt, wp storeFactor, std::vector<Grid<wp>>& uStore){

};

void comNSSolver::updateVars(Grid<wp>&dt, wp storeFactor, std::vector<Grid<wp>>& uStore, std::vector<Grid<wp>>& rStore){

};

wp comNSSolver::getResNorm(){

};

void comNSSolver::getResidual(){

};

bool comNSSolver::isConverged(){
        return false;
};

void comNSSolver::reconstructStates(){


};

void comNSSolver::computeFluxes(){


};

std::vector<wp> comNSSolver::conservativeToPrimitive(std::vector<wp>qc){
        std::vector<wp>qp;
        qp.resize(4);

        wp cvT;

        qp[0]=qc[0];
        qp[1]=qc[1]/qc[0];
        qp[2]=qc[2]/qc[0];
        cvT=(qc[3]/qc[0])-0.5*(qp[1]*qp[1]+qp[2]*qp[2]);
        qp[3]=qp[0]*cvT*(params["GM"]-1.0);

        return qp;

};

std::vector<wp> comNSSolver::primitiveToConservative(std::vector<wp>qp){
        std::vector<wp>qc;
        qc.resize(4);

        wp T;

        qc[0]=qp[0];
        qc[1]=qp[0]*qp[1];
        qc[2]=qp[0]*qp[2];
        T=qp[3]/(params["RGAS"]*qp[0]);
        qc[3]=qp[0]*((params["RGAS"]/(params["GM"]-1.0))*T+0.5*(qp[1]*qp[1]+qp[2]*qp[2]));

        return qc;

};

