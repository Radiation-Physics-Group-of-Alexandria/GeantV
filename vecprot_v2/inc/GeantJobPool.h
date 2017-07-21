#ifndef GEANTV_GEANTJOBPOOL_H
#define GEANTV_GEANTJOBPOOL_H


#include <json.hpp>
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

class GeantHPCJobPool {
public:
  virtual GeantHPCJob GetJob(int n) = 0;
  virtual void ReturnToPool(GeantHPCJob job) = 0;
  virtual bool IsEmpty() = 0;
};

class GeantHepMCJobPool : public GeantHPCJobPool {
public:
  GeantHPCJob GetJob(int n) override;

  void ReturnToPool(GeantHPCJob job) override;

  bool IsEmpty() override;

  void LoadFromFile(std::string fname);

  ~GeantHepMCJobPool(){};
private:
  int JobCounter;
  std::vector<GeantHepMCJob> pool;
};
}}

#endif //GEANTV_GEANTJOBPOOL_H
