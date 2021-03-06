#ifndef __TSPAS_SolidRotationTestCase__
#define __TSPAS_SolidRotationTestCase__

#include "AdvectionTestCase.h"

namespace tspas {

class SolidRotationTestCase : public AdvectionTestCase {
protected:
    double angleSpeed, U0, alpha;
    SpaceCoord *axisPole, *c0, *cr0;
    double R, H0;
public:
    SolidRotationTestCase();
    virtual ~SolidRotationTestCase();

    virtual void init(const ConfigManager &configManager,
                      TimeManager &timeManager);

    Time getStartTime() const;
    Time getEndTime() const;
    double getStepSize() const;

    virtual void calcInitCond(AdvectionManager &advectionManager);

    void calcSolution(double dt, const TimeLevelIndex<2> &timeIdx,
                      AdvectionManager &advectionManager);

    virtual void advance(double time, const TimeLevelIndex<2> &timeIdx);
protected:
    void calcSolution(double dt, const TimeLevelIndex<2> &timeIdx,
                      ScalarField &q);
};

}

#endif // __TSPAS_SolidRotationTestCase__
