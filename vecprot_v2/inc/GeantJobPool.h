#ifndef GEANTV_GEANTJOBPOOL_H
#define GEANTV_GEANTJOBPOOL_H


#include <json.hpp>
#include <chrono>
#include <queue>

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


struct GeantHPCJob {
  int uid;
  std::vector<GeantHepMCJob> hepMCJobs;
};

struct GeantHPCWorker{
  int id;
  int assignedJobId;
  GeantHPCJob assignedJob;
  std::string reqSocket;
  std::chrono::time_point<std::chrono::system_clock> lastContact;
};

class GeantHPCJobPool {
public:
  virtual GeantHPCJob GetJob(int n, const GeantHPCWorker& worker) = 0;
  virtual void ReturnToPool(GeantHPCJob job) = 0;
  virtual bool IsEmpty() = 0;
};

class GeantHepMCJobPool : public GeantHPCJobPool {
public:
  GeantHPCJob GetJob(int n, const GeantHPCWorker& worker) override;

  void ReturnToPool(GeantHPCJob job) override;

  bool IsEmpty() override;

  void LoadFromFile(std::string fname);

  ~GeantHepMCJobPool(){};
private:
  int JobCounter;
  std::vector<std::string> filesForWorkers;
  std::map<std::string,std::vector<GeantHepMCJob>> mapPool;
  std::map<int,std::string> workerFiles;
};
}}

#endif //GEANTV_GEANTJOBPOOL_H
