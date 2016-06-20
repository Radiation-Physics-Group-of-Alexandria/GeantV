#include <cassert>

#include "GUFieldPropagatorPool.h"

// For implementation
#include "GUFieldPropagator.h"

#include <iostream>

// static
std::vector<GUFieldPropagator*> // GUFieldPropagation:: // namespace ...
   GUFieldPropagatorPool::fFieldPropagatorVec;

/// --------------  GUFieldPropagatorPool ------------------------------------
// #include "GUFieldPropagatorPool.h"   // For now, not a separate file

// static
GUFieldPropagatorPool* 
GUFieldPropagatorPool::Instance()
{
   // A lock is REQUIRED for the next line - TODO
   static GUFieldPropagatorPool sInstance;

   return &sInstance;
}

GUFieldPropagatorPool::GUFieldPropagatorPool( GUFieldPropagator* prototype )
   : fInitialisedRKIntegration(false),
     fNumberPropagators(0),
     fPrototype(prototype)
{
   // prototype can be null initially
}

GUFieldPropagatorPool::~GUFieldPropagatorPool()
{
   delete fPrototype;
}

bool
GUFieldPropagatorPool::RegisterPrototype( GUFieldPropagator* prototype )
{
   bool ok = ((fNumberPropagators > 0) && (prototype!=fPrototype)); 
   if( !ok){
      std::cerr << "WARNING from GUFieldPropagatorPool:  "
                << "Changing prototype propagator after having created "
                << fNumberPropagators << " instances. " << std::endl;
   }
   assert( prototype );
   fPrototype= prototype;

   fInitialisedRKIntegration=true;
   return ok;
}

bool
GUFieldPropagatorPool::Initialize( unsigned int numThreads )
{
   if( ! fPrototype ){
       std::cerr << "ERROR> from GUFieldPropagatorPool::Initialize:  "
                 << "Must register prototype propagator before calling Initialize. "
                 << std::endl
                 << "    # propagators= " << fNumberPropagators << " instances. "
                 << std::endl;
       exit(1);
   }
   bool goodExpansion= true;
   if( numThreads > fNumberPropagators )
   {
      // std::cout << "GUFieldPropagatorPool::Initialize  calling Extend for "
      //        << numThreads - fNumberPropagators << " new propagators. " << std::endl;
      Extend( numThreads - fNumberPropagators );
   }

   size_t revSize= fFieldPropagatorVec.size();
   std::cout << " Pool:  revised size= " << revSize
             << " request= " << numThreads << std::endl;

   goodExpansion = ( fFieldPropagatorVec.size() >= numThreads );
   assert (goodExpansion);

   fNumberPropagators= revSize;
   
   return (fPrototype != 0) && goodExpansion;
}

void
GUFieldPropagatorPool::Extend(size_t noNeeded)
{
    size_t num= fFieldPropagatorVec.size();
    assert( fPrototype );
    assert( num < noNeeded );

    while ( num++ < noNeeded )
    {
      //  if( (banks != 0) && (banks[num]!=0) )
      //  fFieldPropagatorVec.push( new(banks(num)) GUFieldPropagator() );
      //  else
      // fFieldPropagatorVec.push_back( new GUFieldPropagator() );
      fFieldPropagatorVec.push_back( fPrototype->Clone() );       
    }
}


#if 0

//// ---------------------  Postpone handling of multiple 
GUFieldPropagator* 
GUFieldPropagatorPool::CreateOrFind( int noNeeded ) // , void** banks )
{
  static int numberCreated= -1;
  static GUFieldPropagatorPool* pInstance= Instance();

  // A lock is REQUIRED for this section - TODO
  if( numberCreated < noNeeded)
  {
    Extend(noNeeded);
    assert( fFieldPropagatorVec.size() == noNeeded );
    // fNum = fFieldPropagatorVec.size();
    numberCreated= noNeeded;
  }
}
#endif
