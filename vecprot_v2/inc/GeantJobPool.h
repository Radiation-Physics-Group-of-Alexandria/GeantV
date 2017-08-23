//===--- GeantJobPool.h - Geant-V --------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantJobPool.h
 * @brief Implementation of job pool for HPC prototype
 * @author Vitalii Drohan
 */
//===----------------------------------------------------------------------===//
#ifndef GEANTV_GEANTJOBPOOL_H
#define GEANTV_GEANTJOBPOOL_H


#include <json.hpp>
#include <chrono>
#include <memory>
#include <set>
#include <Geant/Config.h>
#include "zmq_util.h"
using nlohmann::json;

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

struct GeantHepMCJob {
  std::string fFilename;
  int fOffset;
  int fAmount;
};

void to_json(json& j, const GeantHepMCJob& p);
void from_json(const json& j, GeantHepMCJob& p);


enum class JobType{
  HepMC, Generator
};
struct GeantHPCJob {
  int fUID;
  int fWorkerID;
  JobType fType;
  std::vector<GeantHepMCJob> fHepMCJobs; //Used by HEPMC job
  int fEvents = 0; //Used by Generator job
  ZmqTimer fDispatchTime;
  std::shared_ptr<std::set<int>> fDublicateUIDs;

  GeantHPCJob() : fDublicateUIDs{std::make_shared<std::set<int>>()} {}
};

struct GeantHPCWorker{
  int fID;
  std::string fZMQID;
  std::vector<size_t> fPendingRequests;
  int fDiscardedMsg;
  ZmqTimer fLastContact;
  std::chrono::milliseconds fExpectedTimeForEvent = std::chrono::milliseconds(0);
  int fTransportedEvents = 0;
};

class GeantHPCJobPool {
public:
  virtual GeantHPCJob GetJob(int n, const GeantHPCWorker& worker) = 0;
  virtual void ReturnToPool(GeantHPCJob job) = 0;
  virtual bool IsEmpty() = 0;

  GeantHPCJob GetDublicateJob(GeantHPCJob& job);
protected:
  int GetNextId(){ fJobCounter++; return fJobCounter;};
  int fJobCounter;
};

class GeantHepMCJobPool : public GeantHPCJobPool {
public:
  GeantHPCJob GetJob(int n, const GeantHPCWorker& worker) override;

  void ReturnToPool(GeantHPCJob job) override;

  bool IsEmpty() override;

  void LoadFromFile(std::string fname);

  ~GeantHepMCJobPool() = default;
private:
  std::vector<std::string> filesForWorkers;
  std::map<std::string,std::vector<GeantHepMCJob>> mapPool;
  std::map<int,std::string> workerFiles;
};

class GeantGeneratorJobPool : public GeantHPCJobPool {
public:
  GeantHPCJob GetJob(int n, const GeantHPCWorker& worker) override;

  void ReturnToPool(GeantHPCJob job) override;

  bool IsEmpty() override;

  void SetEventAmount(int ev);

  ~GeantGeneratorJobPool() = default;
private:
  int eventsToDispatch;
};
}}

#endif //GEANTV_GEANTJOBPOOL_H
