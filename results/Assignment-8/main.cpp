#include"icoNSSolver.h"
#include"integrator.h"
#include"interface.h"
Logger logger;


int main(void){
   
        icoNSSolver solver;

        Interface inter(solver);
        Integrator I;
        I.maxSteps=5000;
        I.dumpSteps=5000;
        I.addSolver(solver,"NS");

        return 0;
}
