#include"icoNSSolver.h"
#include"integrator.h"
#include<fenv.h>
Logger logger;


int main(void){
    
        icoNSSolver solver;

        Integrator I;
        I.maxSteps=3;
        I.dumpSteps=3;
        I.addSolver(solver,"NS");
        I.integrate();
        return 0;
}
