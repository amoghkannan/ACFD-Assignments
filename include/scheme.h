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
        extrapolate,
        invalidBCType
};

BCType StringToBCType(std::string str);

enum schemeKey{
        timeStepping,
        timeSteppingLocality,
        derivativeX,
        derivativeY,
        derivativeXX,
        derivativeYY,
        derivativeXY,
        multigrid,
        faceReconstruction,
        stateReconstruction,
        useLimiter,
        fluxCalc,
        viscosityLaw,
        PVCoupling,
        interfaceConfig,
        invalidSchemeKey
};

schemeKey StringToschemeKey(std::string str);

enum schemeVal{
        EULER,
        RK4,
        LOCAL,
        GLOBAL,
        LINEARSOLVER,
        FOU, //First-order upwind
        QUICK,
        C4,   //Fourth-order compact,
        V_CYCLE,
        W_CYCLE,
        MUSCL,
        AUSM,
        CONSTANT,
        SUTHERLANDLAW,
        SIMPLE,
        PISO,
        IBM,
        VOF,
        LEVELSET,
        YES,
        NO,
        NONE,
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
template<typename T>
static T averageValue(T varL, T varR, Mesh&mesh, int i, int j, std::pair<int,int>flags){
        wp volL,volR;

        volL=mesh.volumes(i-flags.first,j-flags.second);

        volR=mesh.volumes(i,j);

        T ans;
        ans=varR*volL+varL*volR;
        ans/=(volR+volL);
        return ans;

};

static Vec2 greenGaussCellBased(Grid<wp>&var,Mesh&mesh, int i, int j);
static Vec2 adjustedFaceGradient(Grid<wp>&var, Mesh&mesh, int i, int j, std::pair<int,int>flags);
static void calcViscosityConstant(dictionary<std::string,wp>&params, Grid<wp>&viscosity, Grid<wp>&temperature);
static void calcViscositySutherland(dictionary<std::string,wp>&params, Grid<wp>&viscosity, Grid<wp>&temperature);
std::vector<std::pair<int,wp>>FOUReconstruction(wp& massFlux);
std::vector<std::pair<int,wp>>QUICKReconstruction(wp& massFlux);

void defaultReconstruction(Grid<wp>&var, Grid<wp>&stateL,Grid<wp>&stateR, int i, int j, std::pair<int,int>flags);
void musclReconstruction(Grid<wp>&var, Grid<wp>&stateL,Grid<wp>&stateR, int i, int j, std::pair<int,int>flags,bool lim);
wp korenLimiter(wp rNum, wp rDenom, bool lim);
wp vanalbadaLimiter(wp rNum, wp rDenom, bool lim);

//Flux calc schemes
void ausmFlux(dictionary<std::string,wp>&params, std::pair<bool,bool>zeroFluxBC, Grid<Vec2>&normals, 
dictionary<std::string,Grid<wp>>&flux, dictionary<std::string,Grid<wp>>&stateL, dictionary<std::string,Grid<wp>>&stateR, 
std::pair<int,int>flags);

void viscousFlux(dictionary<std::string,wp>&params, Mesh& mesh,Grid<wp>&mu, 
dictionary<std::string,Grid<wp>>&flux, dictionary<std::string,Grid<wp>>&state,Grid<wp>&temperature, 
std::pair<int,int>flags);

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
