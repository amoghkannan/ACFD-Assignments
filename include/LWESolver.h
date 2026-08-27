#include"solver.h"

class LWESolver:public Solver{

private:

wp c=1.0;

public:

LWESolver(Mesh& mesh);

virtual void initialCondition() override;

virtual void QDot() override;

virtual void computeTimeStep(Grid& dt) override;

virtual void updateVars(Grid& dt) override;

virtual wp getResNorm() override;

virtual bool isConverged() override;


};
