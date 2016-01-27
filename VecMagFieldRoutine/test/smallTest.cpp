//===--- BenchmarkTiming.cpp - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file  BenchmarkTiming.cpp
 * @brief Benchmark of implementations of bi-linear interpolation of CMS field
 * @author Ananya
 */
//===----------------------------------------------------------------------===//

#include <iostream>

#include <string>
#include <vector>
#include <ctime>
#include <cmath> //for sqrt
// #include <stdlib.h>
#include <cstdlib>

#include <numeric>
#include <string>
#include <functional>

#include <Vc/Vc>

#include "backend/vc/Backend.h"
// #include "backend/vcfloat/Backend.h"
#include "VcFloatBackend.h"
#include "base/Vector.h"

#include "base/Vector3D.h"
#include "base/SOA3D.h"
#include "base/Global.h"

// #include "MagField.h"

#include "CMSmagField.h"

using namespace std;

typedef vecgeom::Vector3D<float> ThreeVector; //normal Vector3D
typedef vecgeom::Vector3D<vecgeom::kVcFloat::precision_v> ThreeVecSimd_t;
typedef vecgeom::Vector<float> VcVectorFloat;




int main(){

    CMSmagField m1;

    m1.ReadVectorData("../VecMagFieldRoutine/cms2015.txt");
    vecgeom::Vector3D<float> position, xyzField;

    std::cout<<"Give x,y and z"<<std::endl;
    std::cin>>position.x()>>position.y()>>position.z();

    m1.GetFieldValue<vecgeom::kScalarFloat>(position, xyzField);

    std::cout<<"Magnetic Field at "<<position<<" is: "<<xyzField<<" tesla."<<std::endl;

}


