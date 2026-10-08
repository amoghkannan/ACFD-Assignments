#include"vars.h"
#include"mesh.h"
#include"LinearSolvers/Jacobi.h"
#include"LinearSolvers/GaussSeidel.h"
#include"LinearSolvers/SGS.h"
#include"LinearSolvers/ILU.h"
#include"LinearSolvers/Cholesky.h"
#include"LinearSolvers/GMRES.h"
#include"LinearSolvers/steepestDescent.h"
#include"LinearSolvers/conjugateGradient.h"
#include"LinearSolvers/biCGStab.h"
#include<fstream>

Logger logger;

int imx=101;
int jmx=101;
int bufW=1;
int bufE=1;
int bufS=1;
int bufN=1;

Grid<wp> phi(Box(imx,jmx,bufW,bufE,bufS,bufN),{NODE,NODE});

std::vector<boundMatEntry> getA(int i,int j,int iLim, int jLim){
        int imx=iLim;
        int jmx=jLim;

        std::vector<boundMatEntry> ans;
        boundMatEntry temp;

        if(i>1 && i<imx && j>1 && j<jmx){
                temp.first={i-1,j};
                temp.second=-1.0;
                ans.push_back(temp);
        
                temp.first={i+1,j};
                temp.second=-1.0;
                ans.push_back(temp);
        
                temp.first={i,j-1};
                temp.second=-1.0;
                ans.push_back(temp);
        
                temp.first={i,j+1};
                temp.second=-1.0;
                ans.push_back(temp); 
        
                temp.first={i,j};
                temp.second=4.0;
                ans.push_back(temp);

        }
        else{
                temp.first={i,j};
                temp.second=1.0;
                ans.push_back(temp);
        };


        return ans;
};

wp getRHS(int i, int j,int iLim, int jLim){
        int imx=iLim;
        int jmx=jLim;

        if(i>1 && i<imx && j>1 && j<jmx){
                return -2.0;
        }
        else{
                return 1.0;
        };


};

int main(void){
        
        wp L=1.0;
        wp delta=L/(imx-1);

        Mesh mesh(Box(imx,jmx,bufW,bufE,bufS,bufN),{NODE,NODE});

        for(int j=0;j<=jmx+1;j++){
                for(int i=0;i<=imx+1;i++){
                        mesh(i,j).x = (i-1)*delta;
                        mesh(i,j).y = (j-1)*delta;
                };
        };


        phi.initVal(0.0);
        std::ofstream outfile("result.dat");

        biCGStab solver(10000,1E-15,1.5,phi,imx,jmx);
        solver.setMatFunc(getA);
        solver.setRHSFunc(getRHS);
//        solver.setupPreconditioner(CHOLESKY_SOLVER,SPLIT_PRECONDITIONER);
        solver.solve(phi);
        phi.print(outfile);
        outfile.close();
        return 0;
}
