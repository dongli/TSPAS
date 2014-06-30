#ifndef __TSPAS_AdvectionManager__
#define __TSPAS_AdvectionManager__

#include "tspas_commons.h"

namespace tspas {

class AdvectionManager {
protected:
    const Domain *domain;
    const Mesh *mesh;
    IOManager io;
    int outputFileIdx;
    
    vector<ScalarField*> Q;
    SingleScalarField Qstar;
    SingleScalarField FX, FY, A, B;
    SingleScalarField volume;

    vec dlon, dlat;
    vec cosLatFull, cosLatHalf;
public:
    AdvectionManager();
    ~AdvectionManager();

    void init(const Domain &domain, const Mesh &mesh,
              const ConfigManager &configManager,
              const TimeManager &timeManager);

    void registerTracer(const string &name, const string &units,
                        const string &brief);

    vector<ScalarField*>& getDensities() { return Q; }

    void input(const TimeLevelIndex<2> &timeIdx, double *q);

    void output(const TimeLevelIndex<2> &newTimeIdx, int ncId);

    void diagnose(const TimeLevelIndex<2> &timeIdx);

    void advance(double dt, const TimeLevelIndex<2> &newTimeIdx,
                 const VelocityField &velocity);
};

}

#endif // __TSPAS_AdvectionManager__
