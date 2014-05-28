#include "tspas.h"

int main(int argc, const char *argv[])
{
    // -------------------------------------------------------------------------
    // check arguments
    if (argc != 2) {
        REPORT_ERROR("Wrong usage! A configuration file path is needed.");
    }
    // -------------------------------------------------------------------------
    geomtk::ConfigManager configManager;
    tspas::AdvectionTestCase *testCase;
    tspas::AdvectionManager advectionManager;
    geomtk::TimeManager timeManager;
    geomtk::TimeLevelIndex<2> oldTimeIdx;
    // -------------------------------------------------------------------------
    // parse configuration
    configManager.parse(argv[1]);
    // -------------------------------------------------------------------------
    // choose test case
    bool isTrueSolution = false;
    std::string testCaseName, subcaseName = "";
    configManager.getValue("test_case", "case_name", testCaseName);
    if (testCaseName == "rotation") {
        testCase = new tspas::SolidRotationTestCase();
        if (configManager.hasKey("test_case", "is_true_solution")) {
            configManager.getValue("test_case", "is_true_solution", isTrueSolution);
        }
    } else if (testCaseName == "deform") {
        testCase = new tspas::DeformationTestCase();
        if (configManager.hasKey("test_case", "subcase")) {
            configManager.getValue("test_case", "subcase", subcaseName);
            testCase->selectSubcase(subcaseName);
        }
    } else if (testCaseName == "barotropic") {
        testCase = new tspas::BarotropicTestCase();
    } else {
        REPORT_ERROR("Unknown test_case \"" << testCaseName << "\"!");
    }
    // -------------------------------------------------------------------------
    // initialization
    timeManager.init(testCase->getStartTime(), testCase->getEndTime(),
                     testCase->getStepSize());
    testCase->init(configManager, timeManager);
    advectionManager.init(testCase->getDomain(), testCase->getMesh(),
                          configManager, timeManager);
    
    testCase->calcInitCond(advectionManager);
    advectionManager.output(oldTimeIdx);
    testCase->advance(timeManager.getSeconds(), oldTimeIdx);
    // -------------------------------------------------------------------------
    // integration loop
    while (!timeManager.isFinished()) {
        geomtk::TimeLevelIndex<2> newTimeIdx = oldTimeIdx+1;
        double time = timeManager.getSeconds()+timeManager.getStepSize();
        testCase->advance(time, newTimeIdx);
        advectionManager.advance(timeManager.getStepSize(), newTimeIdx,
                                 testCase->getVelocityField());
        if (isTrueSolution) {
            testCase->calcSolution(timeManager.getStepSize(),
                                   newTimeIdx, advectionManager);
        }
        timeManager.advance();
        oldTimeIdx.shift();
        advectionManager.output(oldTimeIdx);
    }
    delete testCase;
    return 0;
}