//===--- NumaBlock.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file NumaBlock.h
 * @brief NUMA-aware data blocks
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//
#ifndef GEANT_NUMA_BLOCK
#define GEANT_NUMA_BLOCK

#include <atomic>
#include <cassert>
#include <iostream>
#include <type_traits>
#include "Geant/Config.h"
#include "GeantNuma.h"

namespace Geant {
/**
 * @brief Class NumaBlock
 * @detailed A NUMA block is a concurrent templated factory holding contiguously a
 *           number of objects. It keeps a counters for numbers of used objects.
 */
template <typename T, bool D=false> class NumaBlock {

  using size_t = std::size_t;
  using atomic_size_t = std::atomic<std::size_t>;
  using NumaBlock_t = NumaBlock<T,D>;
  static size_t const cacheline_size = 64;
  typedef char cacheline_pad_t[cacheline_size];

private:
  atomic_size_t fCurrent;  // Current free track index
  atomic_size_t fUsed;     // Number of tracks used

  cacheline_pad_t pad0_;   //! Padding to protect the other data from the hot cache line above

  size_t        fSize;     // Number of tracks stored by the block
  int           fNode;     // NUMA node id
  int           fMaxdepth; // Max depth in case D is used
  NumaBlock_t  *fAddress;  // Address of the block retreivable relative to the one of fArray
  T            *fArray;    //! Array of elements
  

private:
  /** @brief Constructor */
  NumaBlock(size_t size, int node) : fCurrent(0), fUsed(0), fSize(size), fNode(node), fMaxdepth(0), fAddress(this)
  {
    // NUMA block constructor. If the system is NUMA-aware, the block will be alocated
    // on the memory associated with the given NUMA node.
    static_assert(!std::is_polymorphic<T>::value, "Cannot use polymorphic types as NumaBlock");
    static_assert(std::is_default_constructible<T>::value, "Type used in NumaBlock must have default ctor.");
    static_assert(std::is_copy_constructible<T>::value, "Type used in NumaBlock must be copy constructible");
    fArray = reinterpret_cast<T*>(&fArray + 1);
  }

  /** @brief Constructor providing a maximum depth parameter*/
  NumaBlock(size_t size, int node, int maxdepth) : fCurrent(0), fUsed(0), fSize(size), fNode(node), fMaxdepth(maxdepth), fAddress(this)
  {
    // NUMA block constructor. If the system is NUMA-aware, the block will be alocated
    // on the memory associated with the given NUMA node.
    static_assert(!std::is_polymorphic<T>::value, "Cannot use polymorphic types as NumaBlock");
    static_assert(std::is_default_constructible<T>::value, "Type used in NumaBlock must have default ctor.");
    static_assert(std::is_copy_constructible<T>::value, "Type used in NumaBlock must be copy constructible");
    fArray = reinterpret_cast<T*>(&fArray + 1);
    auto el_size = T::SizeOfInstance(maxdepth);
    for (size_t i=0; i<size; ++i) T::MakeInstanceAt((char*)&fArray[0] + i*el_size, maxdepth);
//    std::cout << "Created block: " << this << std::endl;
  }
  
  NumaBlock(const NumaBlock&) = delete;
  NumaBlock& operator=(const NumaBlock&) = delete;

public:
  static NumaBlock *MakeInstance(size_t nvalues, int numa_node)
  {
    // Make an instance. To be released using ReleaseInstance. 
    size_t needed = SizeOfInstance(nvalues);
    void *ptr = numa_aligned_malloc(needed, numa_node, 64);
    NumaBlock *block = new (ptr) NumaBlock(nvalues, numa_node);
    return ( block );    
  }

  static NumaBlock *MakeInstance(size_t nvalues, int numa_node, int maxdepth)
  {
    // Make an instance. To be released using ReleaseInstance. 
    size_t needed = SizeOfInstance(nvalues, maxdepth);
    void *ptr = numa_aligned_malloc(needed, numa_node, 64);
    NumaBlock *block = new (ptr) NumaBlock(nvalues, numa_node, maxdepth);
    return ( block );    
  }

  static void ReleaseInstance(NumaBlock *block)
  {
    // Release the instance of the block
//    std::cout << "deleting block: " << block << std::endl;
    numa_aligned_free(block);
  }
    
  static constexpr size_t SizeOfInstance(size_t nvalues)
  { return ( sizeof(NumaBlock_t) + nvalues*sizeof(T) ); }

  static constexpr size_t SizeOfInstance(size_t nvalues, int maxdepth)
  { return ( sizeof(NumaBlock_t) + nvalues * T::SizeOfInstance(maxdepth) ); }
  
public:
  GEANT_INLINE T *GetValues() { return &fArray[0]; }
  GEANT_INLINE const T *GetValues() const { return &fArray[0]; }

  GEANT_INLINE T &operator[](size_t index) {
    if (!D) return GetValues()[index];
    auto el_size = T::SizeOfInstance(fMaxdepth);
    return *(T*)((char*)&fArray[0]+index*el_size);
  };

  GEANT_INLINE const T &operator[](size_t index) const {
    if (!D) return GetValues()[index];
    auto el_size = T::SizeOfInstance(fMaxdepth);
    return *(T*)((char*)&fArray[0]+index*el_size);
  };
  

  /** @brief Destructor */
  ~NumaBlock() = delete;
  
  /** @brief Get an object pointer from the container */
  GEANT_INLINE T *GetObject(size_t &index) {
    index = fCurrent.fetch_add(1);
    if (index >= fSize) return nullptr;
    fUsed++;
    return ( &operator[](index) );
  }
  
  /** @brief Release an object to the container */
  GEANT_INLINE bool ReleaseObject() { 
    if (fCurrent.load() < fSize) {
      fUsed--;
      return false;
    }
    if ( fUsed.fetch_sub(1) == 1 ) return true;
    return false;
  }

  /** @brief Release an object to the container */
  GEANT_INLINE void Clear() { fCurrent.store(0); fUsed.store(0);}

   /** @brief Returns number of contained objects */
  GEANT_INLINE size_t size() const { return fSize; }

   /** @brief Getter for NUMA node */
  GEANT_INLINE int GetNode() { return fNode; }
 
  /** @brief Check if the block is still in use */
  GEANT_INLINE bool IsDistributed() const { return (fCurrent.load() >= fSize); }
};

} // Geant

#endif