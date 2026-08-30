#pragma once
#include"utils.h"
#include"mesh.h"
#include"vars.h"
#include<string.h>
#include<vector>

enum BCType{
        periodic,
        dirichlet,
        neumann
};

enum schemeKey{
        timeStepping,
        derivativeX,
        derivativeY,
        derivativeXX,
        derivativeYY,
        derivativeXY,
};

enum schemeVal{
        EULER,
        RK4,
        FOU, //First-order upwind
        C4,   //Fourth-order compact
        INVALID
};


class Scheme{

private:

schemeVal invalidScheme=INVALID;
std::vector<schemeKey> keys;
std::vector<schemeVal> values;

BCType BC[4];
wp BCVal[4];

//TDMA coefficients
wp *a=nullptr;
wp *b=nullptr;
wp *c=nullptr;
wp *RHS=nullptr;
wp *x=nullptr;
wp *y=nullptr;

public:

schemeVal getScheme(schemeKey key);
void setScheme(schemeKey key, schemeVal val);
void setBC(int ind, BCType type, wp val);
std::pair<BCType,wp>getBC(int ind);

//Derivative formulas

void TDMA(int n, wp* x);

void EBD1(Grid& phi, Grid& derivative, Mesh& mesh, char dirFlag); //Explicit backward difference (1st order)
void EFD1(Grid& phi, Grid& derivative, Mesh& mesh, char dirFlag); //Explicit forward difference (1st order)
void ECD2(Grid& phi, Grid& derivative, Mesh& mesh, char dirFlag); // Explicit central difference (2nd order)
void ECD2NonUniform(Grid& phi, Grid& derivative, Mesh& mesh, char dirFlag); //2nd order CD for non-uniform grid

void fillCoeffs(Grid& phi, char dirFlag); //Prepare TDMA
void ICD4(Grid& phi, Grid& derivative, Mesh& mesh, char dirFlag); //Implicit (compact) central difference (4th order)

~Scheme(){
        if(a!=nullptr) delete[] a;
        if(b!=nullptr) delete[] b;
        if(c!=nullptr) delete[] c;
        if(RHS!=nullptr) delete[] RHS;
        if(x!=nullptr) delete[] x;
        if(y!=nullptr) delete[] y;

        logger.log("Debug: destroying scheme",1);
};

};
