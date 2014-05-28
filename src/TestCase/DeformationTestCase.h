#ifndef __TSPAS_DeformationTestCase__
#define __TSPAS_DeformationTestCase__

#include "AdvectionTestCase.h"

namespace tspas {

class DeformationTestCase : public AdvectionTestCase {
protected:
    double period;
public:
    DeformationTestCase();
    ~DeformationTestCase();

    virtual void init(const ConfigManager &configManager,
                      const TimeManager &timeManager);

    Time getStartTime() const;
    Time getEndTime() const;
    double getStepSize() const;

    void calcInitCond(AdvectionManager &advectionManager);
    void advance(double time, const TimeLevelIndex<2> &timeIdx);
};

}

#endif // __TSPAS_DeformationTestCase__