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

Grid<wp> Solver::restrictVar(Mesh& fineMesh, Grid<wp>& fineVar){
        Box coarseBox;
        coarseBox.imx=(fineVar.box.imx+1)/2;
        coarseBox.jmx=(fineVar.box.jmx+1)/2;
        coarseBox.bufW=fineVar.box.bufW;
        coarseBox.bufE=fineVar.box.bufE;
        coarseBox.bufS=fineVar.box.bufS;
        coarseBox.bufN=fineVar.box.bufN;

        if(fineVar.indexType[0]==NODE && fineVar.indexType[1]==NODE){

                Grid<wp> coarseVar(coarseBox,{NODE,NODE});
        
                for(int j=1;j<=fineVar.box.jmx;j+=2){
                        for(int i=1;i<=fineVar.box.imx;i+=2){
                                coarseVar((i+1)/2,(j+1)/2)=fineVar(i,j);
                        };
                        
                };
        
                return coarseVar;

        }
        else if(fineVar.indexType[0]==CELL && fineVar.indexType[1]==CELL){
                Grid<wp> coarseVar(coarseBox,{CELL,CELL});
                if(fineMesh.volumes.trueSize==0) logger.log("Warning: Volumes uninitialized",1);
        
                for(int j=1;j<=fineVar.intVect[1];j+=2){
                        for(int i=1;i<=fineVar.intVect[0];i+=2){
                                coarseVar((i+1)/2,(j+1)/2)=volumeAverage(fineMesh,fineVar,i,j);
                        };
                        
                };
        
                return coarseVar;
        }
        else{
                logger.log("Warning: face-based variables not supported for multigrid",1);
        };
};

wp Solver::volumeAverage(Mesh&fineMesh,Grid<wp>&fineVar,int i,int j){

        wp ans=fineVar(i,j)*fineMesh.volumes(i,j);
        wp weights=fineMesh.volumes(i,j);

        ans=ans+fineVar(i+1,j)*fineMesh.volumes(i+1,j);
        weights=weights+fineMesh.volumes(i+1,j);

        ans=ans+fineVar(i,j+1)*fineMesh.volumes(i,j+1);
        weights=weights+fineMesh.volumes(i,j+1);

        ans=ans+fineVar(i+1,j+1)*fineMesh.volumes(i+1,j+1);
        weights=weights+fineMesh.volumes(i+1,j+1);

        return ans/weights;
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
                for(int j=1;j<coarseVar.intVect[1];j++){
                        for(int i=1;i<coarseVar.intVect[0];i++){

                                fineVar(2*i,2*j)+=bilinearInterp4({cc(i,j,coarseMesh),cc(i+1,j,coarseMesh),
                                                                   cc(i+1,j+1,coarseMesh),cc(i,j+1,coarseMesh)},
                                                                  {coarseVar(i,j),coarseVar(i+1,j),
                                                                   coarseVar(i+1,j+1),coarseVar(i,j+1)},
                                                                   cc(2*i,2*j,fineMesh));
                                fineVar(2*i+1,2*j)+=bilinearInterp4({cc(i,j,coarseMesh),cc(i+1,j,coarseMesh),
                                                                   cc(i+1,j+1,coarseMesh),cc(i,j+1,coarseMesh)},
                                                                  {coarseVar(i,j),coarseVar(i+1,j),
                                                                   coarseVar(i+1,j+1),coarseVar(i,j+1)},
                                                                   cc(2*i+1,2*j,fineMesh));
                                fineVar(2*i+1,2*j+1)+=bilinearInterp4({cc(i,j,coarseMesh),cc(i+1,j,coarseMesh),
                                                                   cc(i+1,j+1,coarseMesh),cc(i,j+1,coarseMesh)},
                                                                  {coarseVar(i,j),coarseVar(i+1,j),
                                                                   coarseVar(i+1,j+1),coarseVar(i,j+1)},
                                                                   cc(2*i+1,2*j+1,fineMesh));
                                fineVar(2*i,2*j+1)+=bilinearInterp4({cc(i,j,coarseMesh),cc(i+1,j,coarseMesh),
                                                                   cc(i+1,j+1,coarseMesh),cc(i,j+1,coarseMesh)},
                                                                  {coarseVar(i,j),coarseVar(i+1,j),
                                                                   coarseVar(i+1,j+1),coarseVar(i,j+1)},
                                                                   cc(2*i,2*j+1,fineMesh));
                        };
                };

        };

};

void Solver::restriction(std::vector<Grid<wp>>&R_store){
        for(int i=0;i<nVars;i++){
                R_store[i]=varsDot(i);
        };

        getResNorm();
        
        for(auto key:vars.getKeys()){
                coarser->vars[key]=restrictVar(*mesh,vars[key]);
                coarser->vars[key].initVal(0.0);
                coarser->varsDot[key]=restrictVar(*mesh,varsDot[key]);
        };

        (coarser->ls)->setRHSFunc(coarser->varsDot(0));
};

void Solver::prolongation(std::vector<Grid<wp>>&R_store){
        
        for(auto key:vars.getKeys()){
                prolongateVar(*(finer->mesh),(*mesh),finer->vars[key],vars[key]);
        };

        for(int i=0;i<nVars;i++){
                finer->varsDot(i)=R_store[i];
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

        wp x_xi=nds[1].x-nds[0].x;
        wp y_xi=nds[1].y-nds[0].y;
        wp x_eta=nds[3].x-nds[0].x;
        wp y_eta=nds[3].y-nds[0].y;

        wp delta_x=nd.x-nds[0].x;
        wp delta_y=nd.y-nds[0].y;
        wp J=x_xi*y_eta-y_xi*x_eta;

        wp delta_xi=(y_eta*delta_x-x_eta*delta_y)/J;
        wp delta_eta=(-y_xi*delta_x+x_xi*delta_y)/J;


        return vals[0]*(1.0-delta_xi)*(1.0-delta_eta)+
               vals[1]*(delta_xi)*(1.0-delta_eta)+
               vals[2]*(delta_xi)*(delta_eta)+
               vals[3]*(1.0-delta_xi)*(delta_eta);
};

Node Solver::cc(int i, int j, Mesh& mesh){

        return (mesh(i,j)+mesh(i+1,j)+mesh(i+1,j+1)+mesh(i,j+1))*0.25;
};
