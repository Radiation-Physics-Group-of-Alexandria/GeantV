#include "GeantDistributeManager.h"
#include <mpi.h>
#include "Geant/Error.h"
#include <iostream>

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {
GeantDistributeManger::GeantDistributeManger(GeantRunManager *runManager, GeantConfig *config)
    : fRunMgr(runManager), fConfig(config)
{
}

void GeantDistributeManger::InitializeDistributedApplication(int argc, char *argv[])
{
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    assert("MPI_Init failed");
    std::abort();
  }

  int world_rank;
  int world_size;

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  fConfig->fNClients = world_size;

  hostnames = new char[world_size][256];
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  MPI_Allgather(hostname, 256, MPI_CHAR, hostnames, 256, MPI_CHAR, MPI_COMM_WORLD);

  if (world_rank == 0) {
    fMaster = new GeantEventDispatcher(fConfig);
  }

  fWorker = new GeantEventReceiver(std::string(hostnames[0]), fConfig, fRunMgr);
  fRunMgr->SetEventReceiver(fWorker);

  MPI_Finalize();
}

void GeantDistributeManger::RunDistributedSimulation()
{
  std::thread eventDispThread;
  if (fMaster != nullptr) {
    eventDispThread = std::thread([&] {
      fMaster->Initialize();
      fMaster->RunReqReplyLoop();
    });
  }
  if (fWorker != nullptr) {
    fWorker->Initialize();
    fWorker->Run();
  }
  if (fMaster != nullptr) {
    eventDispThread.join();
  }
}

GeantDistributeManger::~GeantDistributeManger()
{
  delete[] hostnames;
}
}
}
