#include"solver.h"

Solver::Solver(){
        logger.log("Debug: Base solver initialization",1);
};

void Solver::setMesh(Mesh& meshIn){
       mesh=&meshIn;
};

void Solver::setMesh(Box box_){
        mesh=new Mesh(box_, {NODE,NODE});
};

void Solver::setVar(std::string name, std::array<idtype,2>indexType_){
       vars[name]=Grid<wp>(mesh->box,indexType_);
       varsDot[name]=Grid<wp>(mesh->box,indexType_);
       nVars=nVars+1;
};

void Solver::setVar(int ind, Grid<wp>& VarIn){
       vars(ind)=VarIn;
};

void Solver::setScheme(Scheme& schemeIn){
        scheme=schemeIn;
};

void Solver::setScheme(schemeKey key, schemeVal val){
        scheme.setScheme(key,val);
};

void Solver::setBC(std::string boundary, BCType type, wp val){
        scheme.setBC(boundary,type,val);
};

Mesh* Solver::getMesh(){
        return mesh;
};

Scheme& Solver::getScheme(){
        return scheme;
};

void Solver::restrictVar(Mesh& fineMesh, Mesh& coarseMesh, Grid<wp>& fineVar, Grid<wp>&coarseVar){

        if(fineVar.indexType[0]==NODE && fineVar.indexType[1]==NODE){

                for(int j=1;j<=coarseVar.intVect[1];j++){
                        for(int i=1;i<=coarseVar.intVect[0];i++){
                                coarseVar(i,j)=fineVar(2*i-1,2*j-1);
                        };
                        
                };
        
        }
        else if(fineVar.indexType[0]==CELL && fineVar.indexType[1]==CELL){
        
                for(int j=1;j<=coarseVar.intVect[1];j++){
                        for(int i=1;i<=coarseVar.intVect[0];i++){
                                coarseVar(i,j)=volumeAverage(fineMesh,fineVar,i,j)/coarseMesh.volumes(i,j);
                        };
                };
        
        }
        else{
                logger.log("Warning: face-based variables not supported for multigrid",1);
        };

};

wp Solver::volumeAverage(Mesh&fineMesh,Grid<wp>&fineVar,int iCoarse,int jCoarse){

        int i=iCoarse;
        int j=jCoarse;

        wp ans=fineVar(2*i-1,2*j-1)*fineMesh.volumes(2*i-1,2*j-1);

        ans=ans+fineVar(2*i,2*j-1)*fineMesh.volumes(2*i,2*j-1);

        ans=ans+fineVar(2*i-1,2*j)*fineMesh.volumes(2*i-1,2*j);

        ans=ans+fineVar(2*i,2*j)*fineMesh.volumes(2*i,2*j);

        return ans;
};

void Solver::prolongateVar(Mesh& fineMesh, Mesh& coarseMesh, Grid<wp>&fineVar, Grid<wp>& coarseVar){

        if(fineVar.indexType[0]==NODE && fineVar.indexType[1]==NODE){
                for(int j=1;j<=coarseVar.intVect[1];j++){
                        for(int i=1;i<=coarseVar.intVect[0];i++){

                                fineVar(2*i-1,2*j-1)+=coarseVar(i,j);
                                if(2*i-1>=1 && 2*j-2>=1){
                                        fineVar(2*i-1,2*j-2)+=bilinearInterp2({fineMesh(2*i-1,2*j-3),fineMesh(2*i-1,2*j-1)},
                                                                   {coarseVar(i,j-1),coarseVar(i,j)},fineMesh(2*i-1,2*j-2));
                                };

                                if(2*i-2>=1 && 2*j-1>=1){
                                        fineVar(2*i-2,2*j-1)+=bilinearInterp2({fineMesh(2*i-3,2*j-1),fineMesh(2*i-1,2*j-1)},
                                                                   {coarseVar(i-1,j),coarseVar(i,j)},fineMesh(2*i-2,2*j-1));
                                };

                                if(2*i-2>=1 && 2*j-2>=1){
                                        fineVar(2*i-2,2*j-2)+=bilinearInterp4({fineMesh(2*i-3,2*j-3),fineMesh(2*i-1,2*j-3),
                                                                  fineMesh(2*i-1,2*j-1),fineMesh(2*i-3,2*j-1)},
                                                                  {coarseVar(i-1,j-1),coarseVar(i,j-1),
                                                                   coarseVar(i,j),coarseVar(i-1,j)},fineMesh(2*i-2,2*j-2));
                                };
                        };
                };

        }
        
        else if(fineVar.indexType[0]==CELL && fineVar.indexType[1]==CELL){
                for(int j=1;j<=fineVar.intVect[1];j++){
                        for(int i=1;i<=fineVar.intVect[0];i++){
                               
                                int iCoarse=i/2;
                                int jCoarse=j/2;

                                Node nd1=coarseMesh.cc(iCoarse,jCoarse);
                                wp val1=coarseVar(iCoarse,jCoarse);
                                Node nd2=coarseMesh.cc(iCoarse+1,jCoarse);
                                wp val2=coarseVar(iCoarse+1,jCoarse);
                                Node nd3=coarseMesh.cc(iCoarse+1,jCoarse+1);
                                wp val3=coarseVar(iCoarse+1,jCoarse+1);
                                Node nd4=coarseMesh.cc(iCoarse,jCoarse+1);
                                wp val4=coarseVar(iCoarse,jCoarse+1);

                                fineVar(i,j)+=bilinearInterp4({nd1,nd2,nd3,nd4},{val1,val2,val3,val4},fineMesh.cc(i,j));
                       
                        };
                };

        };

};

void Solver::restriction(){

        if(getScheme().getScheme(multigrid)==FAS){

                for(auto key:vars.getKeys()){
                        restrictVar(*mesh,*(coarser->mesh),vars[key],coarser->vars[key]);
                };

                coarser->QDot();

                int counter=0;
                for(auto key:vars.getKeys()){
                        Grid<wp>&coarseVarsStore=coarser->res_store[counter];
                        Grid<wp>&coarseVarsDot=coarser->varsDot[key];
                        Grid<wp>&fineVarsDot=varsDot[key];
                        for(int j=1;j<=coarseVarsDot.intVect[1];j++){
                                for(int i=1;i<=coarseVarsDot.intVect[0];i++){
                                        coarseVarsStore(i,j)=fineVarsDot(2*i-1,2*j-1)+fineVarsDot(2*i-1,2*j)
                                                           +fineVarsDot(2*i,2*j-1)+fineVarsDot(2*i,2*j);
                                        coarseVarsStore(i,j)-=coarseVarsDot(i,j);
                                };
                        };
                        counter++;
                };

        }
        else{
                res_store.clear(); //TODO: Seems inefficient to clear every time step?

                for(int i=0;i<nVars;i++){
                        res_store.push_back(varsDot(i));
                };
                
                getResidual();
                
                for(auto key:vars.getKeys()){
                        coarser->vars[key].initVal(0.0);
                        restrictVar(*mesh,*(coarser->mesh),varsDot[key],coarser->varsDot[key]);
                };

                (coarser->ls)->setRHSFunc(coarser->varsDot(0));
       };

};

void Solver::prolongation(){

        for(int n=0;n<nVars;n++){
                Grid<wp>&currVar=vars(n);
                Grid<wp>&currUStore=u_store[n];
                for(int j=0;j<=currVar.intVect[1]+1;j++){
                        for(int i=0;i<=currVar.intVect[0]+1;i++){
                                if(i==0 || i==currVar.intVect[0]+1 || j==0 || j==currVar.intVect[1]+1){
                                        currUStore(i,j)=0.0;
                                }
                                else{
                                        currUStore(i,j)=currVar(i,j)-currUStore(i,j);
                                };
                        };
                };

                prolongateVar(*(finer->mesh),(*mesh),finer->vars(n),currUStore);
        
        };

        if(getScheme().getScheme(multigrid)!=FAS){
                for(int i=0;i<nVars;i++){
                        finer->varsDot(i)=finer->res_store[i];
                };
        };

};

wp Solver::bilinearInterp2(std::array<Node,2>nds,std::array<wp,2>vals,Node nd){
        wp wt0,wt1,wtS;
        
        Vec2 v0(nd,nds[0]);
        Vec2 v1(nd,nds[1]);
        
        wt0=v0.norm2();
        wt1=v1.norm2();
        wtS=wt0+wt1;

        wt0=wt0/wtS;
        wt1=wt1/wtS;

        return wt0*vals[1]+wt1*vals[0];

};

wp Solver::bilinearInterp4(std::array<Node,4>nds,std::array<wp,4>vals,Node nd){

        //Refer https://iquilezles.org/articles/ibilinear/

        nds[0]=nds[0]-nds[0];
        nds[1]=nds[1]-nds[0];
        nds[2]=nds[2]-nds[0];
        nds[3]=nds[3]-nds[0];

        wp alpha=nds[2].x-nds[1].x-nds[3].x;
        wp beta=nds[2].y-nds[1].y-nds[3].y;

        wp coeffsA=alpha*nds[3].y-beta*nds[2].x;
        wp coeffsB=nd.x*beta-nd.y*alpha-nds[3].x*nds[1].y+nds[1].x*nds[3].y;
        wp coeffsC=nd.x*nds[1].y-nd.y*nds[1].x;

        wp xi,eta;
        wp eta1,eta2;

        wp D;

        if(fabs(coeffsA)<=10.0*SMALL){
                eta=-coeffsC/coeffsB;
        }
        else{
                D=coeffsB*coeffsB-4.0*coeffsA*coeffsC;
                if(D>=0.0){
                        eta1=0.5*(-coeffsB+sqrt(D))/(coeffsA+sign(coeffsA)*SMALL);
                        eta2=0.5*(-coeffsB-sqrt(D))/(coeffsA+sign(coeffsA)*SMALL);
                        eta=(eta1>=0.0 && eta1<=1.0)?eta1:eta2;
                }
                else{
                        return 0.25*(vals[0]+vals[1]+vals[2]+vals[3]);
                };
        };      

        wp num,denom;
        num=nd.y-eta*nds[3].y;
        denom=nds[1].y+eta*(nds[2].y-nds[1].y-nds[3].y);
        
        if(fabs(denom)<10.0*SMALL){
                num=nd.x-eta*nds[3].x;
                denom=nds[1].x+eta*(nds[2].x-nds[1].x-nds[3].x);
        };

        xi=num/denom;

        return  vals[0]*(1.0-xi)*(1.0-eta)
                +vals[1]*(xi)*(1.0-eta)
                +vals[2]*(xi)*(eta)
                +vals[3]*(1.0-xi)*(eta);

};


