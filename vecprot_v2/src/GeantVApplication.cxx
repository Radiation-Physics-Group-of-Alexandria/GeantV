#include "GeantVApplication.h"
#ifdef USE_ROOT
ClassImp(GeantVApplication)
#endif
    //______________________________________________________________________________
    GeantVApplication::GeantVApplication()
#ifdef USE_ROOT
    : TObject() 
#endif
    {
  // Ctor..

}
