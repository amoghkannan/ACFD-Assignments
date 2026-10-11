#include"icoNSSolver.h"
#include"integrator.h"
#include<fenv.h>
Logger logger;


int main(void){
    
        icoNSSolver solver;

        Integrator I;
        I.maxSteps=5000;
        I.dumpSteps=5000;
        I.addSolver(solver,"NS");
        I.integrate();
        return 0;
}
