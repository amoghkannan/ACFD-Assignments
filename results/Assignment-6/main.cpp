#include"icoNSSolver.h"
#include"integrator.h"

Logger logger;


int main(void){
     
        icoNSSolver solver;

        Integrator I;
        I.maxSteps=1;
        I.dumpSteps=1;
        I.addSolver(solver,"NS");
        I.integrate();
        return 0;
}
