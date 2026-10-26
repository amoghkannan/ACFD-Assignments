#include"comNSSolver.h"
#include"integrator.h"

Logger logger;

int main(void){
   
        comNSSolver solver;

        Integrator I;
        I.maxSteps=1;
        I.dumpSteps=500;
        I.addSolver(solver,"NS");
        I.integrate();

        return 0;
}
