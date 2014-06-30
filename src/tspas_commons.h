#ifndef __tspas_commons__
#define __tspas_commons__

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <armadillo>
#include "geomtk.h"

namespace tspas {

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::string;
using std::vector;
    
using arma::vec;

const int FULL = geomtk::RLLStagger::GridType::FULL;
const int HALF = geomtk::RLLStagger::GridType::HALF;
const int CENTER = geomtk::RLLStagger::Location::CENTER;
const int X_FACE = geomtk::RLLStagger::Location::X_FACE;
const int Y_FACE = geomtk::RLLStagger::Location::Y_FACE;
const int FULL_DIMENSION = geomtk::RLLSpaceDimensions::FULL_DIMENSION;;
    
using geomtk::PI2;
using geomtk::RAD;

typedef geomtk::SphereDomain Domain;
typedef geomtk::SphereCoord SpaceCoord;
typedef geomtk::RLLMesh Mesh;
typedef geomtk::NumericRLLField<double, 2> ScalarField;
typedef geomtk::NumericRLLField<double, 1> SingleScalarField;
typedef geomtk::TimeManager TimeManager;
typedef geomtk::IOManager<geomtk::RLLDataFile> IOManager;
typedef geomtk::IOFrequencyUnit IOFrequencyUnit;
template <int N>
using TimeLevelIndex = geomtk::TimeLevelIndex<N>;
typedef geomtk::Time Time;
typedef geomtk::TimeUnit TimeUnit;
typedef geomtk::ConfigManager ConfigManager;
typedef geomtk::RLLVelocityField VelocityField;

}

#endif // __tspas_commons__
