#include"interface.h"

Interface::Interface(Mesh& meshIn){
        readInInterface();
        readMesh(meshIn);

        phi=Grid<wp>(mesh->box,{CELL,CELL});

        computeLevelSet();

        std::ofstream outfile("levelSet.dat");
        for(int j=1;j<=phi.intVect[1];j++){
                for(int i=1;i<=phi.intVect[0];i++){
                        if(phi(i,j)==BIG) continue;
                        outfile<<(*mesh)(i,j).x<<"\t"<<(*mesh)(i,j).y<<"\t"<<phi(i,j)<<std::endl;
                };
                outfile<<std::endl;
        };

        outfile.close();
};

void Interface::readInInterface(){
        std::ifstream infile(ibFile);
        infile>>nIB;

        wp x1,y1,x2,y2,nx,ny;

        iblines.resize(nIB);

        for(int i=0;i<nIB;i++){

                infile>>x1>>y1>>x2>>y2>>nx>>ny;
                iblines[i]=interfaceLine(Vec2(Node(x1,y1),Node(x2,y2)),Vec2(nx,ny));
        };

        infile.close();

};

void Interface::readMesh(Mesh& mesh_){
        mesh=&mesh_;
};

IBDistStat Interface::IBLineDist(Node nd, int idx){

        IBDistStat ans;

        wp alpha;
        wp beta;
        wp touchX,touchY;
        wp dist;

        wp det11,det12,det21,det22;
        det11=iblines[idx].l.en.x-iblines[idx].l.st.x;
        det12=-iblines[idx].n.en.x;
        det21=iblines[idx].l.en.y-iblines[idx].l.st.y;
        det22=-iblines[idx].n.en.y;

        wp b1,b2;
        b1=nd.x-iblines[idx].l.st.x;
        b2=nd.y-iblines[idx].l.st.y;

        wp det=det11*det22-det12*det21;

        alpha=(det22*b1-det12*b2)/det;
        beta=(-det21*b1+det11*b2)/det;
        touchX=iblines[idx].l.st.x+alpha*(iblines[idx].l.en.x-iblines[idx].l.st.x);
        touchY=iblines[idx].l.st.y+alpha*(iblines[idx].l.en.y-iblines[idx].l.st.y);
        dist=Vec2(nd,Node(touchX,touchY)).norm2();
        
        ans=IBDistStat(alpha,beta,dist,Node(touchX,touchY));
        return ans;
};

wp Interface::nearestIBLineDist(Node nd){

        IBDistStat ans;
        
        IBDistStat currAns(0.0,0.0,BIG,Node(0.0,0.0));
        int currIdx=-1;

        for(int i=0;i<=nIB;i++){
                ans=IBLineDist(nd,i);
                if(ans.dist<currAns.dist){
                        currIdx=i;
                        currAns=ans;
                };
        };

        if(currAns.alpha<=0.0 || currAns.alpha>=1.0){ //Need for pseudonormal
                Vec2 nl,v;

                int idx1,idx2;

                idx1=currIdx;

                if(currAns.alpha<=0.0) idx2=(idx1-1+nIB)%nIB;
                if(currAns.alpha>=1.0) idx2=(idx1+1)%nIB;

                nl=(iblines[idx1].n+iblines[idx2].n)*0.5;
                nl/=nl.norm2();

                currAns.alpha=std::min((wp)1.0,std::max(currAns.alpha,(wp)0.0));
                
                if(currAns.alpha==0.0){
                        currAns.touchPoint=iblines[idx1].l.st;
                }
                else if(currAns.alpha==1.0){
                        currAns.touchPoint=iblines[idx1].l.en;
                };

                currAns.dist=Vec2(nd,currAns.touchPoint).norm2();
                
                currAns.beta=nl.dotProduct(Vec2(currAns.touchPoint,nd));
                currAns.beta=currAns.dist*currAns.beta/fabs(currAns.beta);

        };

        return currAns.beta;

};

void Interface::computeLevelSet(){

        wp xMin,xMax,yMin,yMax;

        xMin=BIG;
        xMax=-BIG;
        yMin=BIG;
        yMax=-BIG;

        for(int i=0;i<nIB;i++){
                if(iblines[i].l.st.x<xMin) xMin=iblines[i].l.st.x;
                if(iblines[i].l.st.x>xMax) xMax=iblines[i].l.st.x;
                if(iblines[i].l.st.y<yMin) yMin=iblines[i].l.st.y;
                if(iblines[i].l.st.y>yMax) yMax=iblines[i].l.st.y;
        };

       xMin=0.8*xMin;
       yMin=0.8*yMin;
       xMax=1.2*xMax;
       yMax=1.2*yMax;
       Node cc; 

        for(int j=1;j<=phi.intVect[1];j++){
                for(int i=1;i<=phi.intVect[0];i++){
                       cc=mesh->cc(i,j);
                       if(cc.x>=xMin && cc.x<=xMax && cc.y>=yMin && cc.y<=yMax){
                                phi(i,j)=nearestIBLineDist(cc);
                       }
                       else{
                                phi(i,j)=BIG;
                       };
                };
        };

};

void Interface::initializeVolumeFractionField(){
//Tentative f values(0/1) already set throughout flow field, need to initialize exact values 
//for cells cut by the boundary

};


