
#include "GeantCudaUtils.h"
#include "backend/cuda/Interface.h"

#include "CoprocessorBrokerKernel.h"
#include "GeantTaskData.h"
#include "GeantTrack.h"

namespace Geant {

inline namespace cuda {
template void MakeInstanceArrayAt(GeantTaskData *addr, size_t nElements, size_t sizeOf, size_t, int, GeantPropagator *);

template void MakeInstanceAt(GeantTrack_v *addr, unsigned int, int);

__global__ void Clear(GeantTrack_v *tracks) { tracks->Clear(); }

int Clear_gpu(vecgeom::cxx::DevicePtr<Geant::cuda::GeantTrack_v> &tracks, int blocksPerGrid, int threadsPerBlock,
              cudaStream_t stream)
{
  Clear<<<blocksPerGrid, threadsPerBlock, 0, stream>>>(tracks);
  GEANT_CUDA_ERROR(cudaGetLastError());
  return 1;
}

} // cuda
} // Geant

namespace vecgeom {
namespace cxx {
template void DevicePtr<Geant::cuda::GeantConfig>::Construct() const;
template size_t DevicePtr<Geant::cuda::GeantConfig>::SizeOf();
template void DevicePtr<Geant::cuda::GeantPropagator>::Construct(int) const;
template size_t DevicePtr<Geant::cuda::GeantPropagator>::SizeOf();
template size_t DevicePtr<Geant::cuda::GeantTaskData>::SizeOf();
template size_t DevicePtr<Geant::cuda::GeantTrack_v>::SizeOf();
} // cxx
} // vecgeom
