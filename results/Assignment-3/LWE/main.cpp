#include"integrator.h"
#include"LWESolver.h"

Logger logger;

int main(void){
        
        int imx=100;
        int jmx=1;
        int bufW=1;
        int bufE=1;
        int bufS=0;
        int bufN=0;

        wp L=1.0;
        wp delta=L/(imx-1);

        Mesh mesh(Box(imx,jmx,bufW,bufE,bufS,bufN),{NODE,NODE});

        for(int i=0;i<=imx+1;i++){
                mesh(i,1).x = (i-1)*delta;
                mesh(i,1).y = 0.0;
        };

        LWESolver solver1(mesh);
        solver1.setScheme(timeStepping,RK4);
        solver1.setScheme(derivativeX,FOU);
        solver1.setBC("W",periodic,1);
        solver1.setBC("E",periodic,0);

        LWESolver solver2(mesh);
        solver2.setScheme(timeStepping,RK4);
        solver2.setScheme(derivativeX,C4);
        solver2.setBC("W",periodic,1);
        solver2.setBC("E",periodic,0);

        Integrator I;

        I.setCFL(1.0);
        I.maxSteps=10000;
        I.dumpSteps=1;

        I.addSolver(solver1,"FOU");
        I.addSolver(solver2,"Compact");

        I.integrate();

        return 0;
}
