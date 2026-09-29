#include"mesh.h"
#include"vars.h"
#include"scheme.h"

Logger logger;

int main(void){
        
        int imx=20;
        int jmx=1;
        int bufW=1;
        int bufE=1;
        int bufS=0;
        int bufN=0;

        wp L=1.0;
        wp delta=L/(imx-1);

        Mesh mesh(Box(imx,jmx,bufW,bufE,bufS,bufN),{NODE,NODE});

        mesh(0,1).x=-0.01;
        mesh(1,1).x=0.0;
        mesh(0,1).y=0.0;
        mesh(1,1).y=0.0;
        wp alpha=1.5;

        for(int i=2;i<=imx+1;i++){
                delta=mesh(i-1,1).x-mesh(i-2,1).x;
                mesh(i,1).x = mesh(i-1,1).x+alpha*delta;
                mesh(i,1).y = 0.0;
        };

 
//        for(int i=0;i<=imx+1;i++){
//                mesh(i,1).x = 0.01*(i-1);
//                mesh(i,1).y = 0.0;
//        };

        Grid<wp> phi(Box(imx,jmx,bufW,bufE,bufS,bufN),{NODE,NODE});
        wp currX;

        for(int i=0;i<=imx+1;i++){
                currX=mesh(i,1).x;
                phi(i,1)=currX*currX+5.0*currX-10.0;
        };
        
        Grid<wp> derivative1(Box(imx,jmx,bufW,bufE,bufS,bufN),{NODE,NODE});
        Grid<wp> derivative2(Box(imx,jmx,bufW,bufE,bufS,bufN),{NODE,NODE});

        Scheme scheme;

        scheme.ECD2NonUniform(phi,derivative1,mesh,'x');
        scheme.ECD2(phi,derivative2,mesh,'x');

        std::ofstream outfile("results.dat");

        for(int i=1;i<=imx;i++){
                currX=mesh(i,1).x;
                outfile<<currX<<"\t"<<derivative1(i,1)<<"\t"<<derivative2(i,1)<<"\t"<<2.0*currX+5.0<<std::endl;
        };

        outfile.close();
        return 0;
}
