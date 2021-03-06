//  Create and maintain a 'pool' of Field Propagator instances
//  Each thread will use the same instance, indexed by its thread-id (0<tid<=maxTID)
//

//  For Multi-threaded version only -- ie not for CUDA

//  An implementation on GPU will require a different approach, potentially
//   - a dynamically created set of classes GUFieldPropagator/Driver/Stepper ...
//   - a static method for each class (revised GUFieldPropagator/Driver/..)
//        without the need for an object
//   - a mixture of the two approaches.
//  Current status: to be investigated, once classes are more stable.

#ifndef GUFIELD_PROPAGATOR_POOL_H
#define GUFIELD_PROPAGATOR_POOL_H 1

#include <cstddef>
#include <vector>

// namespace GUFieldPropagation {
// inline namespace GUFIELDPROPAGATION_IMPL_NAMESPACE {

class GUFieldPropagator;
// #include "GUFieldPropagator.h"

class GUFieldPropagatorPool
{
  public:
    // Access methods
    // static GUFieldPropagator* CreateOrFind(int numThreads);
      // It can be called from many threads -- same value must be returned
      //  numThreads must be constant between calls

    static GUFieldPropagatorPool* Instance();

    bool   RegisterPrototype( GUFieldPropagator* prototype );
     // prototype for the propagators for each thread
    
    bool Initialize(unsigned int numThreads); 
     // Create new propagators for each thread !
    
    bool CheckIndex(size_t num){
       assert(num< fFieldPropagatorVec.size());
       // ((void)num); // make compiler happy
       return ( num < fFieldPropagatorVec.size() );
    }
    
    GUFieldPropagator* GetPropagator(int num) {
       CheckIndex(num);       
       return fFieldPropagatorVec[num];
    }

  private:

    GUFieldPropagatorPool( GUFieldPropagator* prototype = 0); // , void** banks=0 );  // Ensure one per thread
    ~GUFieldPropagatorPool();

    void Extend(size_t Num);
     // Create additional propagators, so that total is 'Num'
  private:
    // Invariants -- constant during simulation
    unsigned int fNumberPropagators;
    const  GUFieldPropagator* fPrototype; //  Owned
    
    static std::vector<GUFieldPropagator*> fFieldPropagatorVec;
};

// }
// }
#endif
