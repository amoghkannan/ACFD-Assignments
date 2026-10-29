#include"comNSSolver.h"
#include"interface.h"
#include"integrator.h"

Logger logger;

int main(void){
  
        comNSSolver solver;
        //Interface In(solver); 
        //solver.in=&In;
        
        Integrator I;
        I.setCFL(1.0);
        I.maxSteps=20000;
        I.dumpSteps=1000;
        I.addSolver(solver,"NS");
        I.integrate();

        return 0;
}
