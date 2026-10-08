#include"LinearSolvers/LinearSolver.h"

solverType stringToSolverType(std::string str){
        if(str=="JACOBI") return JACOBI_SOLVER;

        if(str=="GAUSSSEIDEL") return GAUSS_SEIDEL_SOLVER;

        if(str=="SGS") return SGS_SOLVER;

        if(str=="ILU") return ILU_SOLVER;

        if(str=="CHOLESKY") return CHOLESKY_SOLVER;

        if(str=="GMRES") return GMRES_SOLVER;

        if(str=="STEEPESTDESCENT") return STEEPEST_DESCENT_SOLVER;

        if(str=="CONJUGATEGRADIENT") return CONJUGATE_GRADIENT_SOLVER;

        if(str=="BICGSTAB") return BICONJUGATE_GRADIENT_SOLVER;

        logger.log("Invalid linear solver string",1);

        return INVALID_SOLVER;
};

