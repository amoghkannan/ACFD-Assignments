#include"PoissonSolver.h"
#include"integrator.h"

Logger logger;

int imx=255;
int jmx=255;

int main(void){
      
        PoissonSolver ps(imx,jmx,1);
        ps.setScheme(multigrid,V_CYCLE);
        
        Integrator I;
        I.maxSteps=1;
        I.dumpSteps=I.maxSteps;
        I.maxMultigridLevel=4;
        I.timeStepsPerLevel=1;
        
        I.addSolver(ps,"PoissonSolver");
        I.integrate();

        return 0;
}
