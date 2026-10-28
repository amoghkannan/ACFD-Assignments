#include"comNSSolver.h"
#include"interface.h"
#include"integrator.h"

Logger logger;

int main(void){
   
        comNSSolver solver;
        Interface In(solver); 
        solver.in=&In;

        solver.reconstructStates();
        In.resetBandInteriorCells(solver.params,solver.vars["rho"],solver.vars["u"],solver.vars["v"],solver.vars["p"],
                                  solver.x_qp_left,solver.y_qp_left);

        std::ofstream outfile("out.dat");
        for(int j=1;j<=solver.vars["u"].intVect[1];j++){
                for(int i=1;i<=solver.vars["u"].intVect[0];i++){
                        outfile<<i<<"\t"<<j<<"\t"<<In.cellClass(i,j)<<std::endl;        
                };
                outfile<<std::endl;
        };

        outfile.close();

        //Integrator I;
        //I.setCFL(1.0);
        //I.maxSteps=1;
        //I.dumpSteps=1000;
        //I.addSolver(solver,"NS");
        //I.integrate();

        return 0;
}
