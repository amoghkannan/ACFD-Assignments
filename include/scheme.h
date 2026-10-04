#pragma once
#include"utils.h"
#include"mesh.h"
#include"vars.h"
#include<string.h>
#include<vector>

enum BCType{
        periodic,
        dirichlet,
        neumann,
        invalidBCType
};

BCType StringToBCType(std::string str);

enum schemeKey{
        timeStepping,
        derivativeX,
        derivativeY,
        derivativeXX,
        derivativeYY,
        derivativeXY,
        multigrid,
        faceReconstruction,
        PVCoupling,
        invalidSchemeKey
};

schemeKey StringToschemeKey(std::string str);

enum schemeVal{
        EULER,
        RK4,
        LINEARSOLVER,
        FOU, //First-order upwind
        C4,   //Fourth-order compact,
        V_CYCLE,
        W_CYCLE,
        SIMPLE,
        INVALIDSCHEMEVAL
};

schemeVal StringToschemeVal(std::string str);

class Scheme{

private:



//TDMA coefficients
wp *a=nullptr;
wp *b=nullptr;
wp *c=nullptr;
wp *RHS=nullptr;
wp *x=nullptr;
wp *y=nullptr;

public:

dictionary<std::string,std::pair<BCType,wp>> BCDict;
dictionary<schemeKey,schemeVal> schemeDict;

schemeVal& getScheme(schemeKey keyIn);
void setScheme(schemeKey key, schemeVal val);
void setBC(std::string boundary, BCType type, wp val);
std::pair<BCType,wp>getBC(std::string boundary);
bool hasScheme(schemeKey key);

//Derivative formulas

void TDMA(int n, wp* x);
void fillCoeffs(Grid<wp>& phi, char dirFlag); //Prepare TDMA

//Finite difference
void EBD1(Grid<wp>& phi, Grid<wp>& derivative, Mesh& mesh, char dirFlag); //Explicit backward difference (1st order)
void EFD1(Grid<wp>& phi, Grid<wp>& derivative, Mesh& mesh, char dirFlag); //Explicit forward difference (1st order)
void ECD2(Grid<wp>& phi, Grid<wp>& derivative, Mesh& mesh, char dirFlag); // Explicit central difference (2nd order)
void ECD2NonUniform(Grid<wp>& phi, Grid<wp>& derivative, Mesh& mesh, char dirFlag); //2nd order CD for non-uniform grid
void ICD4(Grid<wp>& phi, Grid<wp>& derivative, Mesh& mesh, char dirFlag); //Implicit (compact) central difference (4th order)

//Finite volume
wp greenGaussCellBased(Grid<wp>&var,Mesh&mesh,Vec2 dir, int i, int j);
void faceStatesFOU(Grid<boundMatRow>&coeffs, Grid<wp>&vels, Mesh& mesh, char dir);
void faceStatesQUICK(Grid<boundMatRow>&coeffs, Grid<wp>&vels, Mesh&mesh, char dir);

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
