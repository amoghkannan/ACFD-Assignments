#include"PoissonSolver.h"
#include"integrator.h"

Logger logger;

int imx=257;
int jmx=257;

int main(void){
      
        PoissonSolver ps(imx,jmx,1);
        ps.setScheme(multigrid,W_CYCLE);
        
        Integrator I;
        I.maxSteps=5000;
        I.dumpSteps=I.maxSteps;
        I.maxMultigridLevel=5;
        
        I.addSolver(ps,"PoissonSolver");
        I.integrate();

        return 0;
}
