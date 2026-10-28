#include"interface.h"

Interface::Interface(Solver& s_){
        logger.log("Initializing interface",1);

        s=&s_;
        setMesh(*s->mesh);

        std::ifstream infile("interface.dat");
       
        std::string data,data1,data2;
        wp dataNum;

        logger.log("Reading input files",1);
        while(infile>>data){
                if(data.front()=='#') continue;
                
                if(data=="CONFIG"){
                                logger.log("Reading interface configuration",1);
                                while(infile>>data){
                                        if(data.front()=='#') break;
                                        infile>>data1;
                                        setScheme(StringToschemeKey(data),StringToschemeVal(data1));
                                };
                }
        };

        infile.close();

        readInInterface();
        
        setVar("phi",{CELL,CELL});
        
        if(scheme.getScheme(interfaceConfig)==VOF){
                setVar("f",{CELL,CELL});
                normals=Grid<Vec2>(mesh->box,{CELL,CELL});
                alphas=Grid<wp>(mesh->box,{CELL,CELL});
                F=Grid<wp>(mesh->box,{NODE,CELL});
                G=Grid<wp>(mesh->box,{CELL,NODE});
        }
        else if(scheme.getScheme(interfaceConfig)==IBM){
                vars["phi"].initVal(BIG);
                cellClass=Grid<int>(mesh->box,{CELL,CELL});
                cellClass.initVal(0);
                computeLevelSet();
                classifyCells();
                IDBandCells();
        };

        logger.log("Finished interface setup",1);

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
        det12=iblines[idx].n.en.x;
        det21=iblines[idx].l.en.y-iblines[idx].l.st.y;
        det22=iblines[idx].n.en.y;

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

nearestIBStat Interface::nearestIBLineDist(Node nd){

        IBDistStat ans;
        
        IBDistStat currAns(0.0,0.0,BIG,Node(0.0,0.0));
        int currIdx=-1;

        for(int i=0;i<nIB;i++){
                ans=IBLineDist(nd,i);

                if(ans.dist<currAns.dist){
                        currIdx=i;
                        currAns=ans;
                };
        };

        nearestIBStat ret;
        ret.dist=currAns.dist*currAns.beta;
        ret.touchPoint=currAns.touchPoint;
        ret.nearestIBLineInd=currIdx;
        ret.n=iblines[currIdx].n;

        if(currAns.alpha<0.0 || currAns.alpha>1.0){ //Need for pseudonormal

                Vec2 nl,v;

                int idx1,idx2;

                idx1=currIdx;

                if(currAns.alpha<0.0) idx2=(idx1-1+nIB)%nIB;
                if(currAns.alpha>1.0) idx2=(idx1+1)%nIB;

                IBDistStat altAns=IBLineDist(nd,idx2);

                if(altAns.dist==currAns.dist && altAns.alpha>=0.0 && altAns.alpha<=1.0){ //TODO: Check whether bug in FEST3D too
                        currAns=altAns;
                        ret.dist=currAns.dist*currAns.beta;
                        ret.touchPoint=currAns.touchPoint;
                        ret.nearestIBLineInd=idx2;
                        ret.n=iblines[idx2].n;

                        return ret;
                };

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

                ret.dist=currAns.dist*currAns.beta;
                ret.touchPoint=currAns.touchPoint;
                ret.nearestIBLineInd=currIdx;
                ret.n=nl;

        };

        return ret;
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

        for(int j=1;j<=vars["phi"].intVect[1];j++){
                for(int i=1;i<=vars["phi"].intVect[0];i++){
                       cc=mesh->cc(i,j);
                       if(cc.x>=xMin && cc.x<=xMax && cc.y>=yMin && cc.y<=yMax){
                                nearestIBStat ret=nearestIBLineDist(cc);
                                vars["phi"](i,j)=ret.dist;
                       }
                       else{
                                vars["phi"](i,j)=BIG;
                       };
                };
        };

};

void Interface::classifyCells(){
       
        Grid<wp>&phi=vars["phi"];

        for(int j=1;j<=vars["phi"].intVect[1];j++){
                for(int i=1;i<=vars["phi"].intVect[0];i++){
                        if(phi(i,j)>=0.0) cellClass(i,j)=0;
                        if(phi(i,j)<0.0) cellClass(i,j)=-1;
                };
        };

};

void Interface::IDBandCells(){
       
        Node fc;

        for(int j=1;j<=cellClass.intVect[1];j++){
                for(int i=1;i<=cellClass.intVect[0];i++){
                        
                        bool fieldNeighbour=false;

                        if(cellClass(i,j)==-1){

                                if(cellClass(i-1,j)==0){
                                        fieldNeighbour=true;
                                        fc=((*mesh)(i,j)+(*mesh)(i,j+1))*0.5;
                                        nearestIBStat ret=nearestIBLineDist(fc);
                                        bandfaces.emplace_back(interfaceFace(i,j,i-1,j,i,j,fc,'x',ret));
                                };

                                if(cellClass(i+1,j)==0){
                                        fieldNeighbour=true;
                                        fc=((*mesh)(i+1,j)+(*mesh)(i+1,j+1))*0.5;
                                        nearestIBStat ret=nearestIBLineDist(fc);
                                        bandfaces.emplace_back(interfaceFace(i,j,i+1,j,i+1,j,fc,'x',ret));
                                };
        
                                if(cellClass(i,j-1)==0){
                                        fieldNeighbour=true;
                                        fc=((*mesh)(i,j)+(*mesh)(i+1,j))*0.5;
                                        nearestIBStat ret=nearestIBLineDist(fc);
                                        bandfaces.emplace_back(interfaceFace(i,j,i,j-1,i,j,fc,'y',ret));
                                };

                                if(cellClass(i,j+1)==0){
                                        fieldNeighbour=true;
                                        fc=((*mesh)(i,j+1)+(*mesh)(i+1,j+1))*0.5;
                                        nearestIBStat ret=nearestIBLineDist(fc);
                                        bandfaces.emplace_back(interfaceFace(i,j,i,j+1,i,j+1,fc,'y',ret));
                                };
                
                                if(fieldNeighbour){
                                       cellClass(i,j)=-2;
                                       bandcells.emplace_back(interfaceCell(i,j)); 
                                }
                                else{
                                       interiorcells.emplace_back(interfaceCell(i,j)); 
                                };

                        };
                };
        };

};

IBInterpStat Interface::IBInterp(interfaceFace& f, Grid<wp>& u, Grid<wp>& v){
        IBInterpStat ans;

        int ic=f.fieldcellI;
        int jc=f.fieldcellJ;

        int is=ic-1<1?1:ic-1;
        int js=jc-1<1?1:jc-1;

        int ie=ic+1>u.intVect[0]?u.intVect[0]:ic+1;
        int je=jc+1>u.intVect[1]?u.intVect[1]:jc+1;

        ans.interpU=0.0;
        ans.interpV=0.0;
        ans.delf=0.0;
        wp wtSum=0.0;
        wp currWt;
        Node tp=f.stat.touchPoint;
        Vec2 n=f.stat.n;

        for(int j=js;j<=je;j++){
                for(int i=is;i<=ie;i++){
                        if(cellClass(i,j)<0) continue;
                        Node cc=mesh->cc(i,j);
                        Vec2 vec(cc.x-tp.x,cc.y-tp.y);
                        if(vec.dotProduct(n)<0) continue;
                        wp dPar=vec.dotProduct(n);
                        wp dPer=(vec-n*dPar).norm2();
                        currWt=1.0/(dPer+SMALL);
                        ans.interpU=ans.interpU+u(i,j)*currWt;
                        ans.interpV=ans.interpV+v(i,j)*currWt;
                        ans.delf=ans.delf+dPar*currWt;
                        wtSum=wtSum+currWt;
                };
        };

        ans.interpU=ans.interpU/wtSum;
        ans.interpV=ans.interpV/wtSum;
        ans.delf=ans.delf/wtSum;

        return ans;
};

void Interface::setSlipVelocityFace(Grid<wp>&u, Grid<wp>&v,
dictionary<std::string,Grid<wp>>&x_qp_l,dictionary<std::string,Grid<wp>>&x_qp_r,
dictionary<std::string,Grid<wp>>&y_qp_l,dictionary<std::string,Grid<wp>>&y_qp_r)
{
        Grid<wp>&uxL=x_qp_l["u"];
        Grid<wp>&vxL=x_qp_l["v"];
        Grid<wp>&uxR=x_qp_r["u"];
        Grid<wp>&vxR=x_qp_r["v"];

        Grid<wp>&uyL=y_qp_l["u"];
        Grid<wp>&vyL=y_qp_l["v"];
        Grid<wp>&uyR=y_qp_r["u"];
        Grid<wp>&vyR=y_qp_r["v"];

        for(auto& face:bandfaces){
               int iface=face.bandfaceI;
               int jface=face.bandfaceJ;
               IBInterpStat stat=IBInterp(face,u,v);

               wp uInterp=stat.interpU;
               wp vInterp=stat.interpV;
               wp delf=stat.delf;
               wp dface=face.stat.dist;

               Vec2 velInterp(uInterp,vInterp);
               Vec2 velFace=velInterp-face.stat.n*(velInterp.dotProduct(face.stat.n)*(1.0-(dface/delf)));

               if(face.faceType=='x'){
                        uxL(iface,jface)=velFace.en.x;
                        vxL(iface,jface)=velFace.en.y;
                        uxR(iface,jface)=velFace.en.x;
                        vxR(iface,jface)=velFace.en.y;
               }
               else{
                        uyL(iface,jface)=velFace.en.x;
                        vyL(iface,jface)=velFace.en.y;
                        uyR(iface,jface)=velFace.en.x;
                        vyR(iface,jface)=velFace.en.y;
               };

        };

};

void Interface::setPressureDensityFace(Grid<wp>&rho, Grid<wp>&p,
dictionary<std::string,Grid<wp>>&x_qp_l,dictionary<std::string,Grid<wp>>&x_qp_r,
dictionary<std::string,Grid<wp>>&y_qp_l,dictionary<std::string,Grid<wp>>&y_qp_r)
{
        Grid<wp>&pxL=x_qp_l["p"];
        Grid<wp>&rhoxL=x_qp_l["rho"];
        Grid<wp>&pxR=x_qp_r["p"];
        Grid<wp>&rhoxR=x_qp_r["rho"];

        Grid<wp>&pyL=y_qp_l["p"];
        Grid<wp>&rhoyL=y_qp_l["rho"];
        Grid<wp>&pyR=y_qp_r["p"];
        Grid<wp>&rhoyR=y_qp_r["rho"];

        for(auto& face:bandfaces){
               int iface=face.bandfaceI;
               int jface=face.bandfaceJ;
               int ifield=face.fieldcellI;
               int jfield=face.fieldcellJ;

               if(face.faceType=='x'){
                        pxL(iface,jface)=p(ifield,jfield);
                        rhoxL(iface,jface)=rho(ifield,jfield);
                        pxR(iface,jface)=p(ifield,jfield);
                        rhoxR(iface,jface)=rho(ifield,jfield);
               }
               else{
                        pyL(iface,jface)=p(ifield,jfield);
                        rhoyL(iface,jface)=rho(ifield,jfield);
                        pyR(iface,jface)=p(ifield,jfield);
                        rhoyR(iface,jface)=rho(ifield,jfield);
               };

        };

};

void Interface::resetBandInteriorCells(dictionary<std::string,wp>&params,Grid<wp>&rho,Grid<wp>&u,Grid<wp>&v,Grid<wp>&p,
dictionary<std::string,Grid<wp>>&x_qp_l,dictionary<std::string,Grid<wp>>&y_qp_l)
{

        Grid<wp>&uxL=x_qp_l["u"];
        Grid<wp>&vxL=x_qp_l["v"];

        Grid<wp>&uyL=y_qp_l["u"];
        Grid<wp>&vyL=y_qp_l["v"];

        for(auto& cell:bandcells){
               int icell=cell.i;
               int jcell=cell.j;

               int counter=0;
               
               wp uAvg,vAvg,pAvg,rhoAvg;
               
               uAvg=0.0;
               vAvg=0.0;
               pAvg=0.0;
               rhoAvg=0.0;

               if(cellClass(icell-1,jcell)==0){
                        uAvg=uAvg+2.0*uxL(icell,jcell)-u(icell-1,jcell);
                        vAvg=vAvg+2.0*vxL(icell,jcell)-v(icell-1,jcell);
                        pAvg=pAvg+p(icell-1,jcell);
                        rhoAvg=rhoAvg+rho(icell-1,jcell);
                        counter=counter+1;
               };

               if(cellClass(icell+1,jcell)==0){
                        uAvg=uAvg+2.0*uxL(icell+1,jcell)-u(icell+1,jcell);
                        vAvg=vAvg+2.0*vxL(icell+1,jcell)-v(icell+1,jcell);
                        pAvg=pAvg+p(icell+1,jcell);
                        rhoAvg=rhoAvg+rho(icell+1,jcell);
                        counter=counter+1;
               };

               if(cellClass(icell,jcell-1)==0){
                        uAvg=uAvg+2.0*uyL(icell,jcell)-u(icell,jcell-1);
                        vAvg=vAvg+2.0*vyL(icell,jcell)-v(icell,jcell-1);
                        pAvg=pAvg+p(icell,jcell-1);
                        rhoAvg=rhoAvg+rho(icell,jcell-1);
                        counter=counter+1;
               };

               if(cellClass(icell,jcell+1)==0){
                        uAvg=uAvg+2.0*uyL(icell,jcell+1)-u(icell,jcell+1);
                        vAvg=vAvg+2.0*vyL(icell,jcell+1)-v(icell,jcell+1);
                        pAvg=pAvg+p(icell,jcell+1);
                        rhoAvg=rhoAvg+rho(icell,jcell+1);
                        counter=counter+1;
               };

               u(icell,jcell)=uAvg/counter;
               v(icell,jcell)=vAvg/counter;
               rho(icell,jcell)=rhoAvg/counter;
               p(icell,jcell)=pAvg/counter;

        };

        wp rhoinf=params["RHOINF"];
        wp pinf=params["PINF"];

        for(auto& cell:interiorcells){
               int icell=cell.i;
               int jcell=cell.j;

               u(icell,jcell)=0.0;
               v(icell,jcell)=0.0;
               rho(icell,jcell)=rhoinf;
               p(icell,jcell)=pinf;

        };

};

Node Interface::InfiniteEdgeIntersection(Vec2 l1, Vec2 l2, bool& intersects){
        Node ans;

        wp Det11,Det12,Det21,Det22;

        Det11=l2.en.x-l2.st.x;
        Det12=-(l1.en.x-l1.st.x);
        Det21=l2.en.y-l2.st.y;
        Det22=-(l1.en.y-l1.st.y);

        wp b1,b2;
        b1=l1.st.x-l2.st.x;
        b2=l1.st.y-l2.st.y;

        wp Det=Det11*Det22-Det12*Det21;

        if(fabs(Det)<SMALL){//l1 || l2; TODO: Correct this
                intersects=false;
                return ans;
        };

        intersects=true;

        wp alpha=(Det22*b1-Det12*b2)/Det;

        ans.x=l2.st.x+alpha*(l2.en.x-l2.st.x);
        ans.y=l2.st.y+alpha*(l2.en.y-l2.st.y);

        return ans;
};

wp Interface::shoelaceArea(std::vector<Node>nds){
        int n=nds.size();
        if(n<3) return 0.0;

        nds.push_back(nds[0]);

        wp ans=0.0;

        for(int i=0;i<n;i++){
                ans=ans+nds[i].x*nds[i+1].y;
                ans=ans-nds[i].y*nds[i+1].x;
        };

        return 0.5*fabs(ans);
};

bool Interface::isLeft(Node nd, Vec2 seg){
        Vec2 v(seg.en,nd);

        wp det11=seg.en.x-seg.st.x;
        wp det12=seg.en.y-seg.st.y;
        wp det21=v.en.x-v.st.x;
        wp det22=v.en.y-v.st.y;

        wp checker=det11*det22-det12*det21;

        return checker>0.0;

};

std::vector<Node> Interface::SutherlandHodgman(std::vector<Vec2>&clipPolygon,std::vector<Node>&subjectPolygon){
        std::vector<Node>outputList;
        std::vector<Node>inputList;

        outputList=subjectPolygon;

        bool intersects;

        for(auto edge:clipPolygon){
                inputList=outputList;
                outputList.clear();

                int n=inputList.size();

                for(int i=0;i<n;i++){
                        Node currPoint=inputList[i];
                        Node prevPoint=inputList[(i-1+n)%n];

                        Node intersectPoint=InfiniteEdgeIntersection(Vec2(prevPoint,currPoint),edge,intersects);

                        if(isLeft(currPoint,edge)){
                                if(!isLeft(prevPoint,edge)){
                                        outputList.push_back(intersectPoint);
                                };

                                outputList.push_back(currPoint);
                        }
                        else if(isLeft(prevPoint,edge)){
                                outputList.push_back(intersectPoint);
                        };

                };

        };

        return outputList;
};

void Interface::initializeVolumeFractionField(){

        std::vector<Node>subjectPolygon;
        std::vector<Vec2>clipPolygon;
        std::vector<Node>outputList;

        for(int n=0;n<nIB;n++){
                subjectPolygon.push_back(iblines[n].l.st);
        };

        for(int j=1;j<=vars["f"].intVect[1];j++){
                for(int i=1;i<=vars["f"].intVect[0];i++){
                        if(vars["phi"](i,j)==BIG) continue;

                        clipPolygon.clear();
                        clipPolygon.push_back(Vec2((*mesh)(i,j),(*mesh)(i+1,j)));
                        clipPolygon.push_back(Vec2((*mesh)(i+1,j),(*mesh)(i+1,j+1)));
                        clipPolygon.push_back(Vec2((*mesh)(i+1,j+1),(*mesh)(i,j+1)));
                        clipPolygon.push_back(Vec2((*mesh)(i,j+1),(*mesh)(i,j)));

                        outputList=SutherlandHodgman(clipPolygon,subjectPolygon);
                        wp fCurr=shoelaceArea(outputList)/mesh->volumes(i,j);

                        vars["f"](i,j)=fCurr==0.0?(vars["phi"](i,j)>0.0?0.0:1.0):fCurr;

                };
        };
        
};

Vec2 Interface::calcYoungNormal(int i, int j){
        return Scheme::greenGaussCellBased(vars["f"],*mesh,i,j)*(-1.0);
};

Vec2 Interface::calcHtFuncNormal(int i, int j, char dir){

        wp CL=0.0;
        wp CR=0.0;
        wp hL,hR;
        wp sgn;
        wp m;
        Vec2 ans;

        if(dir=='x'){
                for(int j1=j-1;j1<=j+1;j1++){
                        hL=0.5*(mesh->normalsJ(i-1,j1).norm2()+mesh->normalsJ(i-1,j1+1).norm2());
                        CL+=vars["f"](i-1,j1)*mesh->volumes(i-1,j1)/hL;
                        hR=0.5*(mesh->normalsJ(i+1,j1).norm2()+mesh->normalsJ(i+1,j1+1).norm2());
                        CR+=vars["f"](i+1,j1)*mesh->volumes(i+1,j1)/hR;
                };
                
                sgn=vars["f"](i,j)-vars["f"](i,j+1);
                sgn=sgn>0.0?1.0:-1.0;
                m=-(CR-CL)/(mesh->cc(i+1,j).x-mesh->cc(i-1,j).x);
                m/=sgn;
                return Vec2(m,sgn);
        }
        else if(dir=='y'){
                for(int i1=i-1;i1<=i+1;i1++){
                        hL=0.5*(mesh->normalsI(i1,j-1).norm2()+mesh->normalsI(i1+1,j-1).norm2());
                        CL+=vars["f"](i1,j-1)*mesh->volumes(i1,j-1)/hL;
                        hR=0.5*(mesh->normalsI(i1,j+1).norm2()+mesh->normalsI(i1+1,j+1).norm2());
                        CR+=vars["f"](i1,j+1)*mesh->volumes(i1,j+1)/hR;
                };
                sgn=vars["f"](i,j)-vars["f"](i+1,j);
                sgn=sgn>0.0?1.0:-1.0;
                m=-(CR-CL)/(mesh->cc(i,j+1).y-mesh->cc(i,j-1).y);
                m/=sgn;
                return Vec2(sgn,m);
        };

};

Vec2 Interface::calcNormal(int i, int j){
        Vec2 mYoung=calcYoungNormal(i,j);
        
        Vec2 mxh=calcHtFuncNormal(i,j,'x');
        Vec2 myh=calcHtFuncNormal(i,j,'y');

        if(fabs(mxh.en.x)>fabs(myh.en.y)){ //Compare b/w myh, mYoung
                if(fabs(myh.en.y)>fabs(mYoung.en.y/mYoung.en.x)){
                        return myh;
                }
                else{
                        return mYoung;
                };
        }
        else{ //Compare b/w mxh, mYoung
                if(fabs(mxh.en.y)>fabs(mYoung.en.y/mYoung.en.x)){
                        return mxh;
                }
                else{
                        return mYoung;
                };

        };

};

wp Interface::calcAlpha(int i, int j, Vec2 m){
        wp fCurr=vars["f"](i,j);
        
        bool didSwitch=false;

        if(fCurr>0.5){
                m*=(-1.0);
                fCurr=1.0-fCurr;
                didSwitch=true;
        };

        wp cos_theta=1.0;
        wp sin_theta=0.0;

        if(m.en.x<0.0){
                if(m.en.y<0.0){
                        cos_theta=-1.0;
                        sin_theta=0.0;
                }
                else{
                        cos_theta=0.0;
                        sin_theta=1.0;
                };
        }
        else{
                if(m.en.y<0.0){
                        cos_theta=0.0;
                        sin_theta=-1.0;
                };
        };

        wp temp1,temp2;

        temp1=m.en.x*cos_theta+m.en.y*sin_theta;
        temp2=-m.en.x*sin_theta+m.en.y*cos_theta;

        m.en.x=temp1;
        m.en.y=temp2;
        
        std::array<Node,4>coords;
        coords[0]=(*mesh)(i,j);
        coords[1]=(*mesh)(i+1,j);
        coords[2]=(*mesh)(i+1,j+1);
        coords[3]=(*mesh)(i,j+1);
        
        for(auto& c:coords){
                temp1=c.x*cos_theta+c.y*sin_theta;
                temp2=-c.x*sin_theta+c.y*cos_theta;
                c.x=temp1;
                c.y=temp2;
        };

        wp transX,transY,minX,minY,maxX,maxY;

        minX=std::min({coords[0].x,coords[1].x,coords[2].x,coords[3].x});
        minY=std::min({coords[0].y,coords[1].y,coords[2].y,coords[3].y});

        maxX=std::max({coords[0].x,coords[1].x,coords[2].x,coords[3].x});
        maxY=std::max({coords[0].y,coords[1].y,coords[2].y,coords[3].y});

        transX=minX;
        transY=minY;
        
        for(auto& c:coords){
                c.x-=transX;
                c.y-=transY;
        };

        wp scaleX,scaleY;

        scaleX=fabs(maxX-minX); //TODO: What about curvilinear grids??
        scaleY=fabs(maxY-minY);

        for(auto& c:coords){
                c.x/=scaleX;
                c.y/=scaleY;
        };

        m.en.x/=scaleY;
        m.en.y/=scaleX;

        wp mSum=m.en.x+m.en.y;

        m.en.x/=mSum;
        m.en.y/=mSum;

        wp m2=std::max(fabs(m.en.x),fabs(m.en.y));
        wp m1=std::min(fabs(m.en.x),fabs(m.en.y));

        wp A1=0.5*m1/m2;

        wp alpha;

        if(fCurr<A1){
                alpha=sqrt(2.0*m2*m2*fCurr);
        }
        else{
                alpha=m2*fCurr+0.5*m1;
        };

        alpha*=mSum;
        m.en.x*=mSum;
        m.en.y*=mSum;

        for(auto& c:coords){
                c.x*=scaleX;
                c.y*=scaleY;
        };

        m.en.x*=scaleY;
        m.en.y*=scaleX;
        alpha*=(scaleX*scaleY);
        
        alpha+=(m.en.x*transX+m.en.y*transY);

        if(didSwitch) alpha*=(-1.0);

        return alpha;
};

void Interface::PLIC(){

        for(int j=2;j<normals.intVect[1];j++){
                for(int i=2;i<normals.intVect[0];i++){
                        if(vars["f"](i,j)>0.0 && vars["f"](i,j)<1.0){
                                normals(i,j)=calcNormal(i,j);
                                alphas(i,j)=calcAlpha(i,j,normals(i,j));
                        };
                };
        };

};

void Interface::computeFluxes(wp deltaT){

        std::array<Node,4>cellCoords;

        wp u,v;

        Vec2 currNormal;
        wp currAlpha;

        F.initVal(0.0);
        G.initVal(0.0);

        for(int j=1;j<=vars["f"].intVect[1];j++){
                for(int i=1;i<=vars["f"].intVect[0];i++){
                        if(vars["f"](i,j)==0.0) continue;

                        currNormal=normals(i,j);
                        currAlpha=alphas(i,j);

                        //Left face
                        if((*IFaceContraVel)(i,j)<0.0){
                                
                                u=fabs((*IFaceContraVel)(i,j));

                                cellCoords={(*mesh)(i,j),
                                            Node((*mesh)(i,j).x+u*deltaT,(*mesh)(i,j).y),
                                            Node((*mesh)(i,j+1).x+u*deltaT,(*mesh)(i,j+1).y),
                                            (*mesh)(i,j+1)};

                                
                                F(i,j)=(enclosedArea(currNormal,currAlpha,cellCoords))/deltaT;

                        };

                        //Right face
                        if((*IFaceContraVel)(i+1,j)>0.0){
                                
                                u=fabs((*IFaceContraVel)(i+1,j));

                                cellCoords={Node((*mesh)(i+1,j).x-u*deltaT,(*mesh)(i+1,j).y),
                                            (*mesh)(i+1,j),
                                            (*mesh)(i+1,j+1),
                                            Node((*mesh)(i+1,j+1).x-u*deltaT,(*mesh)(i+1,j+1).y)};
                                                                                                
                                
                                F(i+1,j)=(enclosedArea(currNormal,currAlpha,cellCoords))/deltaT;

                        };

                        //Bottom face
                        if((*JFaceContraVel)(i,j)<0.0){
                                
                                v=fabs((*JFaceContraVel)(i,j));

                                cellCoords={(*mesh)(i,j),
                                            (*mesh)(i+1,j),
                                            Node((*mesh)(i+1,j+1).x,(*mesh)(i+1,j).y+v*deltaT),
                                            Node((*mesh)(i,j+1).x,(*mesh)(i,j+1).y+v*deltaT)};

                                
                                G(i,j)=(enclosedArea(currNormal,currAlpha,cellCoords))/deltaT;

                        };

                        //Top face
                        if((*JFaceContraVel)(i,j+1)>0.0){
                                
                                v=fabs((*JFaceContraVel)(i,j+1));

                                cellCoords={Node((*mesh)(i,j+1).x,(*mesh)(i,j+1).y-v*deltaT),
                                            Node((*mesh)(i+1,j+1).x,(*mesh)(i+1,j+1).y-v*deltaT),
                                            (*mesh)(i+1,j+1),
                                            (*mesh)(i,j+1)};
                                
                                G(i,j+1)=(enclosedArea(currNormal,currAlpha,cellCoords))/deltaT;

                        };

                };
        };

};

bool Interface::isInside(Vec2 m, wp alpha, Node nd){
        wp test=m.en.x*nd.x+m.en.y*nd.y-alpha;

        if(alpha>0){
                return test*alpha>0.0;
        }
        else{
                return test*alpha<0.0;
        };

};

wp Interface::enclosedArea(Vec2 m, wp alpha, std::array<Node,4>cellCoords){
      
        std::vector<Node>cutCellCoords;

        Node currNode,nextNode;
        Node intersectNode;
        Node p1,p2;
        bool intersects;
       
        Vec2 currNormal;
        wp currAlpha;

        for(int i=0;i<=3;i++){
                currNode=cellCoords[i];
                nextNode=cellCoords[(i+1)%4];

                if(isInside(m,alpha,currNode)){
                        cutCellCoords.push_back(currNode);
                       
                        if(!isInside(m,alpha,currNode)){
                                p1=m.en.x!=0.0?Node((alpha/m.en.x),0.0):Node(1.0,alpha/m.en.y);
                                p2=m.en.y!=0.0?Node(0.0,alpha/m.en.y):Node(alpha/m.en.x,1.0);
                                intersectNode=InfiniteEdgeIntersection(Vec2(p1,p2), Vec2(currNode,nextNode), intersects);
                                cutCellCoords.push_back(intersectNode);
                        };

                }
                else{
                        if(isInside(m,alpha,nextNode)){
                                p1=m.en.x!=0.0?Node((alpha/m.en.x),0.0):Node(1.0,alpha/m.en.y);
                                p2=m.en.y!=0.0?Node(0.0,alpha/m.en.y):Node(alpha/m.en.x,1.0);
                                intersectNode=InfiniteEdgeIntersection(Vec2(p1,p2), Vec2(currNode,nextNode), intersects);
                                cutCellCoords.push_back(intersectNode);
                        };

                };
        
        };

        return fabs(shoelaceArea(cutCellCoords));

};

void Interface::initialCondition(){
        vars["f"].initVal(0.0);

        computeLevelSet();

        initializeVolumeFractionField();

        PLIC();

};


void Interface::applyBC(){

};

void Interface::QDot(){

};

void Interface::computeTimeStep(Grid<wp>& dt){

};

void Interface::updateVars(Grid<wp>& dt, wp storeFactor){

};

void Interface::updateVars(Grid<wp>& dt, wp storeFactor, std::vector<Grid<wp>>& uStore){

};

void Interface::updateVars(Grid<wp>& dt, wp storeFactor, std::vector<Grid<wp>>& uStore, std::vector<Grid<wp>>& rStore){

};

wp Interface::getResNorm(){

};

void Interface::getResidual(){

};

bool Interface::isConverged(){
        return s->isConverged();
};

