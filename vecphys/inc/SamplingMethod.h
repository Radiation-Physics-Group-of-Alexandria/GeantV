#ifndef SamplingMethod_H
#define SamplingMethod_H

namespace vecphys {

enum SamplingMethod {
  kNullMethod = -1,
  kAlias,                //Alias method [0]
  kRejection,            //Geant4 CompositionAndRejection [1]
  kUnpack,               //Shuffling [2]
  kNumberSamplingMethod  //[3]
};

} // end namespace vecphys

#endif
