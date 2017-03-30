//===--- GeantPerfEvaluator.h - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantPerfEvaluator.h
 * @brief Performance evaluation class in Geant-V prototype
 * @author Oksana Shadura
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_PERF_EVALUATOR
#define GEANT_PERF_EVALUATOR

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief Class for GA tunning
 * @details [long description]
 * 
 * @param y0 [description]
 * @param y [description]
 * @param e [description]
 * @return [description]
 */
class GeantPerfEvaluator {

private:
  double fMaxResMemory;  /** Maximum resident memory */
  double fMaxVirtMemory; /** Maximum virtual memory */
  double fRunTime;       /** Run time of job */

public:
  GeantPerfEvaluator() : fMaxResMemory(0), fMaxVirtMemory(0), fRunTime(0){};

  ~GeantPerfEvaluator() {}

  double GetMaxResMemory() { return fMaxResMemory; }

  double GetMaxVirtMemory() { return fMaxVirtMemory; }

  double GetRunTimeMemory() { return fRunTime; }

  void SetMaxResMemory(double maxresmem) { fMaxResMemory = maxresmem; }

  void SetMaxVirtMemory(double maxvmem) { fMaxVirtMemory = maxvmem; }

  void SetRunTimeMemory(double rt) { fRunTime = rt; }

  // More functions to arrive...
};

}}

#endif
