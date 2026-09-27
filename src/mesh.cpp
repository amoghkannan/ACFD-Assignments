#include"mesh.h"

void Mesh::calcVolumes(){
        volumes=Grid<wp>(box,{CELL,CELL});
        
        for(int j=1;j<box.jmx;j++){
                for(int i=1;i<box.imx;i++){
                        volumes(i,j)=calcVolume({(*this)(i,j),(*this)(i+1,j),(*this)(i+1,j+1),(*this)(i,j+1)});
                };
        };
};

wp Mesh::calcVolume(std::vector<Node>nodesIn){
        nodesIn.push_back(nodesIn[0]);
        int N=nodesIn.size();

        double ans=0.0;

        for(int n=0;n<N-1;n++){ 
                ans=ans+nodesIn[n].x*nodesIn[n+1].y;
                ans=ans-nodesIn[n].y*nodesIn[n+1].x;
        };

        ans=ans*0.5;
        return fabs(ans);

};

void Mesh::setGhostNodes(){

        for(int j=1;j<=box.jmx;j++){
               for(int i=0;i>=1-box.bufW;i--){
                        (*this)(i,j)=(*this)(i+1,j)*2.0-(*this)(i+2,j);
               };

               for(int i=box.imx+1;i<=box.imx+box.bufE;i++){
                        (*this)(i,j)=(*this)(i-1,j)*2.0-(*this)(i-2,j);
               };

        };

        for(int i=1-box.bufW;i<=box.imx+box.bufE;i++){
                for(int j=0;j>=1-box.bufS;j--){
                        (*this)(i,j)=(*this)(i,j+1)*2.0-(*this)(i,j+2); 
                };

                for(int j=box.jmx+1;j<=box.jmx+box.bufN;j++){
                        (*this)(i,j)=(*this)(i,j-1)*2.0-(*this)(i,j-2); 
                };

        };

};

Mesh Mesh::coarsen(){
        Box coarseBox;
        coarseBox.imx=(int)(box.imx/2)+1;
        coarseBox.jmx=(int)(box.jmx/2)+1;
        coarseBox.bufW=box.bufW;
        coarseBox.bufE=box.bufE;
        coarseBox.bufS=box.bufS;
        coarseBox.bufN=box.bufN;

        Mesh coarseMesh(box,{NODE,NODE});

        for(int j=1;j<=box.jmx;j+=2){
                for(int i=1;i<=box.imx;i+=2){
                        coarseMesh((i+1)/2,(j+1)/2)=(*this)(i,j);
                };
                
        };

        coarseMesh.setGhostNodes();

        return coarseMesh;
};

Mesh::~Mesh(){
logger.log("Debug: destroying mesh",1);
};
