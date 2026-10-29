#include"comNSSolver.h"
#include"interface.h"
#include"integrator.h"
#include<fenv.h>

Logger logger;

int main(void){
  
        feenableexcept(FE_INVALID | FE_OVERFLOW);

        comNSSolver solver;
        //Interface In(solver); 
        //solver.in=&In;
        
        Integrator I;
        I.setCFL(1.0);
        I.maxSteps=1000;
        I.dumpSteps=1000;
        I.addSolver(solver,"NS");
        I.integrate();

        return 0;
}
