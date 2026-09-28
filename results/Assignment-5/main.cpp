#include"PoissonSolver.h"
#include"integrator.h"

Logger logger;

int imx=121;
int jmx=121;

int main(void){
      
        PoissonSolver ps(imx,jmx,1);
        ps.setScheme(multigrid,V_CYCLE);
        
        Integrator I;
        I.maxSteps=500;
        I.dumpSteps=I.maxSteps;
        I.maxMultigridLevel=4;
        
        I.addSolver(ps,"PoissonSolver");
        I.integrate();

        return 0;
}
