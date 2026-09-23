#include"vars.h"
#include"mesh.h"
#include"LinearSolvers/Jacobi.h"
#include"LinearSolvers/GaussSeidel.h"
#include"LinearSolvers/SGS.h"
#include"LinearSolvers/ILU.h"
#include"LinearSolvers/Cholesky.h"
#include"LinearSolvers/GMRES.h"
#include<fstream>

Logger logger;

int imx=101;
int jmx=101;
int bufW=1;
int bufE=1;
int bufS=1;
int bufN=1;

Grid<wp> phi(imx,jmx,bufW,bufE,bufS,bufN);

std::vector<boundMatEntry> getA(int i,int j,Grid<wp>&var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

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

wp getRHS(int i, int j,Grid<wp>&var){
        std::array<int,6>sizeArr=var.size();
        int imx=sizeArr[0];
        int jmx=sizeArr[1];

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

        Mesh mesh(imx,jmx,bufW,bufE,bufS,bufN);

        for(int j=0;j<=jmx+1;j++){
                for(int i=0;i<=imx+1;i++){
                        mesh(i,j).x = (i-1)*delta;
                        mesh(i,j).y = (j-1)*delta;
                };
        };


        phi.initVal(0.5);
        std::ofstream outfile("result.dat");

        GMRES solver(1000,1E-7,1.5,10,100,phi);
        solver.setMatFunc(getA);
        solver.setRHSFunc(getRHS);
        solver.setupPreconditioner(CHOLESKY_SOLVER,SPLIT_PRECONDITIONER);
        solver.solve(phi);
        phi.print(outfile);
        outfile.close();
        return 0;
}
