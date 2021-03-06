#ifndef __TSPAS_BarotropicTestCase__
#define __TSPAS_BarotropicTestCase__

#include "AdvectionTestCase.h"
#include "barotropic_model.h"

namespace tspas {

class BarotropicTestCase : public AdvectionTestCase {
protected:
    barotropic_model::ToyTestCase testCase;
    barotropic_model::BarotropicModel_A_ImplicitMidpoint model;
public:
    BarotropicTestCase();
    virtual ~BarotropicTestCase();

    virtual void init(const ConfigManager &configManager,
                      TimeManager &timeManager);

    Time getStartTime() const;
    Time getEndTime() const;
    double getStepSize() const;

    virtual const Domain& getDomain() const { return model.getDomain(); }
    virtual const Mesh& getMesh() const { return model.getMesh(); }
    
    virtual void calcInitCond(AdvectionManager &advectionManager);
    
    virtual void output(const TimeLevelIndex<2> &timeIdx,
                        AdvectionManager &advectionManager);
    
    virtual void advance(double time, const TimeLevelIndex<2> &timeIdx);
};

}

#endif // __TSPAS_BarotropicTestCase__
