#ifndef ROOT_TNudyAlias
#define ROOT_TNudyAlias
/**
 * @file TNudyAlias.h
 * @brief Header file for Alias Sampling
*/

class TRandom;

#include "TNudyCore.h"


/**
 * @class TNudyAlias
 * @brief Random number generation from a discrete distribution with Alias method
 * @author Harphool Kumawat
*/
class TNudyAlias {
private:
  /** Length of alias table */
  int fLen;
  // Alias table components
  /** Probability density function */
  double *fP;    //[fLen] 
  /** X random variable */
  double *fX;    //[fLen]
  /** Alias table */
  double *fA;    //[fLen]
  /** Residual values */
  double *fR;    //[fLen]
  /** Pointer to uniform random number generator */
  TRandom *fRnd; 
public:
  TNudyAlias();
  virtual ~TNudyAlias();
  TNudyAlias(double *p, double *x, int len, unsigned int seed = 65539);
  void DumpTable() const;
  double Random() const;

  /** @brief Return the probability density function at bin i */
  double GetP(int i) const {
    if (i < fLen)
      return fP[i];
    else
      return -1;
  };
  /** @brief Return the random variable value at bin i */
  double GetX(int i) const {
    if (i < fLen)
      return fX[i];
    else
      return -1;
  };
  /** @brief Return the number of bins of the proability distribution */
  int GetLen() const { return fLen; }

#ifdef TNUDYALIAS_MULTITHREAD
  /**
   * @class TNudyComStruct
   * @brief Class to communicate between thread handler and object
   */
  class TNudyComStruct {
  public:
    /** Array of Alias Tables */
    TNudyAlias *fAl;
    int fI;
    TNudyComStruct(TNudyAlias *a, int j) {
      fAl = a;
      fI = j;
    }
    virtual ~TNudyComStruct() {}
  };
  static void *ThreadHandle(void *ptr);
  double *fMult; //! When generating using a multithreaded approach
  int fMultLen;  //! Number of random values to be generated using the multi threaded approach
  double *Randoms(int n);
#endif

#ifdef USE_ROOT
  ClassDef(TNudyAlias, 1)
#endif
};

#endif
