#ifndef GEANTV_GEANTJOBPOOL_H
#define GEANTV_GEANTJOBPOOL_H


#include <json.hpp>
#include <chrono>
#include <queue>
#include <Geant/Config.h>
#include "zmq_util.h"
#include <memory>
#include <set>
using nlohmann::json;

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

struct GeantHepMCJob {
  std::string filename;
  int offset;
  int amount;
};

void to_json(json& j, const GeantHepMCJob& p);
void from_json(const json& j, GeantHepMCJob& p);


enum class JobType{
  HepMC, Generator
};
struct GeantHPCJob {
  int uid;
  int workerId;
  JobType type;
  std::vector<GeantHepMCJob> hepMCJobs;
  int events = 0;
  ZmqTimer dispatchTime;
  std::shared_ptr<std::set<int>> dublicateUIDs;

  GeantHPCJob() : dublicateUIDs{std::make_shared<std::set<int>>()} {}
};

struct GeantHPCWorker{
  int id;
  std::string zmqID;
  std::vector<size_t> pendingMsg;
  int discardedMsg;
  ZmqTimer lastContact;
  std::chrono::milliseconds expectedTimeForEvent = std::chrono::milliseconds(0);
  int transportedEvents = 0;
};

class GeantHPCJobPool {
public:
  virtual GeantHPCJob GetJob(int n, const GeantHPCWorker& worker) = 0;
  virtual void ReturnToPool(GeantHPCJob job) = 0;
  virtual bool IsEmpty() = 0;

  GeantHPCJob GetDublicateJob(GeantHPCJob& job);
protected:
  int GetNextId(){ JobCounter++; return JobCounter;};
  int JobCounter;
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
