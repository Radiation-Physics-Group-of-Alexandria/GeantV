#include "GeantJobPool.h"

#include <iostream>
#include <fstream>

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

GeantHPCJob GeantHepMCJobPool::GetJob(int n, const GeantHPCWorker& worker) {
  GeantHPCJob res;
  res.fType = JobType::HepMC;
  res.fUID = GetNextId();
  int addedEvents = 0;
  while(addedEvents != n && !mapPool.empty()){
    if(workerFiles.count(worker.fID)>0){
      std::string lastFile = workerFiles[worker.fID];
      auto& filePool = mapPool[lastFile];
      while(addedEvents != n && filePool.size() != 0) {
        res.fHepMCJobs.push_back(*(filePool.end() - 1));
        filePool.pop_back();
        ++addedEvents;
      }
      if(addedEvents!=n){
        workerFiles.erase(worker.fID);
      }
    } else {
      workerFiles[worker.fID] = filesForWorkers.front();
      std::cout << "worker id: " << worker.fID << " file: " << workerFiles[worker.fID] << '\n';
      std::string f = filesForWorkers[0];
      filesForWorkers.erase(filesForWorkers.begin());
      filesForWorkers.push_back(f);
    }
    for(auto it = mapPool.begin(); it != mapPool.end();){
      if(it->second.size()==0){
        it = mapPool.erase(it);
      } else {
        ++it;
      }
    }
  }
  res.fDublicateUIDs->insert(res.fUID);
  return res;
}

bool GeantHepMCJobPool::IsEmpty() {
  return mapPool.empty();
}

void GeantHepMCJobPool::ReturnToPool(GeantHPCJob job) {
  for( auto& j : job.fHepMCJobs){
    auto& filePool = mapPool[j.fFilename];
    auto comp = [](const GeantHepMCJob& lhs, const GeantHepMCJob& rhs){
                                       return lhs.fOffset > rhs.fOffset;
                                     };
    auto sortedIt = std::upper_bound(filePool.begin(),filePool.end(),j, comp);
    filePool.insert(sortedIt,j);
  }
}

void GeantHepMCJobPool::LoadFromFile(std::string fname) {
  std::ifstream inputFile(fname);
  std::string line;
  while(inputFile >> line){ //line format: pathToFile:eventOffset:eventAmount
    auto col1 = line.find(':',0);
    auto col2 = line.find(':',col1+1);

    auto fFileName = line.substr(0, col1);
    auto fOffset = std::stoi(line.substr(col1+1, col2-col1-1));
    auto fEventAmount = std::stoi(line.substr(col2+1));

    for (int i = 0; i < fEventAmount; ++i) {
      GeantHepMCJob job;
      job.fAmount = 1;
      job.fOffset = fOffset;
      job.fFilename = fFileName;
      mapPool[fFileName].push_back(job);
      ++fOffset;
    }
    std::reverse(mapPool[fFileName].begin(),mapPool[fFileName].end());
    filesForWorkers.push_back(fFileName);
  }
}

GeantHPCJob GeantGeneratorJobPool::GetJob(int n, const GeantHPCWorker &worker) {
  GeantHPCJob job;
  job.fType = JobType::Generator;
  job.fEvents = std::min(n,eventsToDispatch);
  job.fUID = GetNextId();
  eventsToDispatch -= job.fEvents;
  return job;
}

bool GeantGeneratorJobPool::IsEmpty() {
  return eventsToDispatch == 0;
}

void GeantGeneratorJobPool::ReturnToPool(GeantHPCJob job) {
  eventsToDispatch += job.fEvents;
}

void GeantGeneratorJobPool::SetEventAmount(int ev) {
  eventsToDispatch = ev;
}

void to_json(json& j, const GeantHepMCJob& p){
  j = json{{"fname",p.fFilename},{"offset",p.fOffset},{"n",p.fAmount}};
}

void from_json(const json& j, GeantHepMCJob& p){
  p.fFilename = j.at("fname").get<std::string>();
  p.fAmount = j.at("n").get<int>();
  p.fOffset = j.at("offset").get<int>();
}

GeantHPCJob GeantHPCJobPool::GetDublicateJob(GeantHPCJob &job) {
  GeantHPCJob newJob = job;
  newJob.fUID = GetNextId();
  newJob.fDublicateUIDs->insert(newJob.fUID);
  return newJob;
}
}
}