#include"icoNSSolver.h"
#include"integrator.h"

Logger logger;


int main(void){
     
        icoNSSolver solver;

        Integrator I;
        I.maxSteps=100;
        I.dumpSteps=10;
        I.addSolver(solver,"NS");
        I.integrate();
        return 0;
}
