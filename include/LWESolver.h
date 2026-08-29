#include"solver.h"

class LWESolver:public Solver{

private:

wp c=1.0;

public:

LWESolver(Mesh& mesh);

void initialCondition() override;

void applyBC() override;

void QDot() override;

void computeTimeStep(Grid& dt) override;

void updateVars(Grid& dt, wp storeFactor) override;

void updateVars(Grid& dt, wp storeFactor, std::vector<Grid>& uStore) override;

void updateVars(Grid& dt, wp storeFactor, std::vector<Grid>& uStore, std::vector<Grid>& rStore) override;

wp getResNorm() override;

bool isConverged() override;

~LWESolver() override{
        logger.log("Debug: Destroying LWEsolver",1);
};

};
