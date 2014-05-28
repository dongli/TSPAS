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
    
    vector<Field*> Q;
    SingleLevelField Qstar;
    SingleLevelField FX, FY, A, B;
    SingleLevelField volume;

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

    void input(const TimeLevelIndex<2> &timeIdx, vector<Field*> &q);

    void output(const TimeLevelIndex<2> &newTimeIdx);

    void diagnose(const TimeLevelIndex<2> &timeIdx);

    void advance(double dt, const TimeLevelIndex<2> &newTimeIdx,
                 const VelocityField &velocity);
};

}

#endif // __TSPAS_AdvectionManager__
