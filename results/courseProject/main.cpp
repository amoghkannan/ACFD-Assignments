#include"comNSSolver.h"
#include"interface.h"
#include"integrator.h"
#ifdef _OPENMP
#include<omp.h>
#endif

Logger logger;

int main(void){
        comNSSolver solver;
        Interface In(solver); 
        solver.in=&In;
        Integrator I;
        I.setCFL(1.0);
        I.maxSteps=20000;
        I.dumpSteps=10000;
        I.maxMultigridLevel=4;
        I.addSolver(solver,"NS");
        I.integrate();

        return 0;
}
