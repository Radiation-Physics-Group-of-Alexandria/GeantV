/// \file scalarfloat/Backend.h
/// temporary version of a scalar backend with main floating point type = float
/// will be superceded by new Backend mechanism that comes with VecCore

#ifndef VECGEOM_BACKEND_SCALARFLOATBACKEND_H_
#define VECGEOM_BACKEND_SCALARFLOATBACKEND_H_

#include "base/Global.h"

#include <algorithm>
#include <cstring>

namespace vecgeom {
inline namespace VECGEOM_IMPL_NAMESPACE {

struct kScalarFloat {
  typedef int       int_v;
  typedef float     precision_v;
  typedef bool      bool_v;
  typedef Inside_t  inside_v;
  const static bool early_returns = true;
  // alternative typedefs ( might supercede above typedefs )
  typedef int                   Int_t;
  typedef Precision       Double_t;
  typedef bool Bool_t;
  typedef int  Index_t; // the type of indices

  constexpr static precision_v kOne = 1.0;
  constexpr static precision_v kZero = 0.0;
  const static bool_v kTrue = true;
  const static bool_v kFalse = false;

  template <class Backend>
  VECCORE_ATT_HOST_DEVICE
  static VECGEOM_CONSTEXPR_RETURN bool IsEqual() { return false; }

  VECCORE_ATT_HOST_DEVICE
  VECGEOM_FORCE_INLINE
  static Precision Convert(Precision const &input) { return input; }
};

#if 0
// ScalarFloat is an auxiliary backend, supposed to be used in conjunction with another
// one (assummedly Scalar), so we can not use this.
#ifdef VECGEOM_SCALAR_FLOAT
constexpr size_t kVectorSize  = 1;
#define VECGEOM_BACKEND_TYPE         vecgeom::kScalarFloat
#define VECGEOM_BACKEND_PRECISION_FROM_PTR(P) (*(P))
#define VECGEOM_BACKEND_PRECISION_TYPE        Precision
//#define VECGEOM_BACKEND_PRECISION_NOT_SCALAR
#define VECGEOM_BACKEND_BOOL         vecgeom::ScalarBool
#define VECGEOM_BACKEND_INSIDE       vecgeom::kScalarFloat::inside_v
#endif
#endif

//template <>
//VECCORE_ATT_HOST_DEVICE
//inline VECGEOM_CONSTEXPR_RETURN bool kScalar::IsEqual<kScalar>() {
//  return true;
//}

typedef kScalar::int_v    ScalarInt;
typedef kScalar::precision_v ScalarDouble;
typedef kScalar::bool_v   ScalarBool;

//template <typename Type>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//void CondAssign(const bool cond,
//                Type const &thenval, Type const &elseval, Type *const output) {
//  *output = (cond) ? thenval : elseval;
//}

//template <typename Type>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//void MaskedAssign(const bool cond,
//                  Type const &thenval, Type *const output) {
//  *output = (cond) ? thenval : *output;
//}

//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//bool IsFull(bool const &cond){
//    return cond;
//}
//
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//bool Any(bool const &cond) {
//  return cond;
//}
//
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//bool IsEmpty(bool const &cond){
//    return !cond;
//}

//template <typename Type>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//Type Pow(Type const &x, Type arg) {
//   return pow(x,arg);
//}
//
//template <typename Type>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//Type Pow(Type const &x, int arg) {
//   return pow(x,arg);
//}
//
//template <typename Type>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//Type Abs(const Type val) {
//  return fabs(val);
//}
//
//template <typename Type>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//Type Sqrt(const Type val) {
//  return std::sqrt(val);
//}
//
//template <typename Type>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//Type Pow(const Type val1, const Type val2) {
//  return std::pow(val1, val2);
//}
//
//template <typename Type>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//Type Cbrt(const Type val1) {
//  return cbrt(val1);
//}
//
//
//template <typename Type>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//Type ATan2(const Type y, const Type x) {
//  if (x != 0) return  std::atan2(y, x);
//  if (y >  0) return  kPi / 2;
//  if (y <  0) return -kPi / 2;
//  return  0;
//}
//
//template <typename T>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//T Min(T const &val1, T const &val2) {
//#ifndef VECGEOM_NVCC_DEVICE
//  return std::min(val1, val2);
//#else
//  return val1 < val2 ? val1 : val2;
//#endif
//}
//
//template <typename T>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//T Max(T const &val1, T const &val2) {
//#ifndef VECGEOM_NVCC_DEVICE
//  return std::max(val1, val2);
//#else
//  return val1 > val2 ? val1 : val2;
//#endif
//}

VECCORE_ATT_HOST_DEVICE
VECGEOM_FORCE_INLINE
float sin(const float radians) {
  return std::sin(radians);
}

VECCORE_ATT_HOST_DEVICE
VECGEOM_FORCE_INLINE
float cos(const float radians) {
  return std::cos(radians);
}

VECCORE_ATT_HOST_DEVICE
VECGEOM_FORCE_INLINE
float tan(const float radians) {
  return std::tan(radians);
}

VECCORE_ATT_HOST_DEVICE
VECGEOM_FORCE_INLINE
float Floor( float val ){
    return std::floor(val);
}

//#ifndef VECGEOM_USOLIDS
//template <typename Type>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//void swap(Type &a, Type &b) {
//  std::swap(a, b);
//}
//#endif
//
//template <typename Type>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//void copy(Type const *begin, Type const *const end, Type *const target) {
//#ifndef VECGEOM_NVCC_DEVICE
//  std::copy(begin, end, target);
//#else
//  std::memcpy(target, begin, sizeof(Type)*(end-begin));
//#endif
//}
//
//template <typename Type>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//void reverse_copy(Type const *const begin, Type const *end,
//                  Type *const target) {
//#ifndef VECGEOM_NVCC_DEVICE
//  std::reverse_copy(begin, end, target);
//#else
//  while (--end >= begin) *target++ = *end;
//#endif
//}
//
//template <typename Type>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//void reverse(Type *begin, Type *end) {
//#ifndef VECGEOM_NVCC_DEVICE
//  std::reverse(begin, end);
//#else
//  while (begin++ < end--) swap(begin, end);
//#endif
//}
//
//template <typename Type>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//Type* AlignedAllocate(size_t size) {
//#ifndef VECGEOM_NVCC
//  return static_cast<Type*>(_mm_malloc(sizeof(Type)*size, kAlignmentBoundary));
//#else
//  return new Type[size];
//#endif
//}
//
//template <typename Type>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//void AlignedFree(Type *allocated) {
//#ifndef VECGEOM_NVCC
//  _mm_free(allocated);
//#else
//  delete[] allocated;
//#endif
//}
//
//template <typename IteratorType>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//IteratorType min_element(IteratorType first, IteratorType last) {
//  return std::min_element(first, last);
//}
//
//template <typename IteratorType>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//bool all_of(IteratorType first, IteratorType last) {
//  return std::all_of(first, last, [](bool b){return b;});
//#endif
//}
//
//template <typename InputIterator1, typename InputIterator2>
//VECCORE_ATT_HOST_DEVICE
//VECGEOM_FORCE_INLINE
//bool equal(InputIterator1 first, InputIterator1 last, InputIterator2 target) {
//#ifndef VECGEOM_NVCC_DEVICE
//  return std::equal(first, last, target);
//#else
//  while (first != last) {
//    if (*first++ != *target++) return false;
//  }
//  return true;
//#endif
//}

} } // End global namespace

#endif // VECGEOM_BACKEND_SCALARFLOATBACKEND_H_
