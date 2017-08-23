#include "GeantJobPool.h"

#include <iostream>
#include <fstream>

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

GeantHPCJob GeantHepMCJobPool::GetJob(int n, const GeantHPCWorker& worker) {
  GeantHPCJob res;
  res.type = JobType::HepMC;
  res.uid = GetNextId();
  int addedEvents = 0;
  while(addedEvents != n && !mapPool.empty()){
    if(workerFiles.count(worker.id)>0){
      std::string lastFile = workerFiles[worker.id];
      auto& filePool = mapPool[lastFile];
      while(addedEvents != n && filePool.size() != 0) {
        res.hepMCJobs.push_back(*(filePool.end() - 1));
        filePool.pop_back();
        ++addedEvents;
      }
      if(addedEvents!=n){
        workerFiles.erase(worker.id);
      }
    } else {
      workerFiles[worker.id] = filesForWorkers.front();
      std::cout << "worker id: " << worker.id << " file: " << workerFiles[worker.id] << '\n';
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
  return res;
}

bool GeantHepMCJobPool::IsEmpty() {
  return mapPool.empty();
}

void GeantHepMCJobPool::ReturnToPool(GeantHPCJob job) {
  for( auto& j : job.hepMCJobs){
    auto& filePool = mapPool[j.filename];
    auto comp = [](const GeantHepMCJob& lhs, const GeantHepMCJob& rhs){
                                       return lhs.offset > rhs.offset;
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
      job.amount = 1;
      job.offset = fOffset;
      job.filename = fFileName;
      mapPool[fFileName].push_back(job);
      ++fOffset;
    }
    std::reverse(mapPool[fFileName].begin(),mapPool[fFileName].end());
    filesForWorkers.push_back(fFileName);
  }
}

GeantHPCJob GeantGeneratorJobPool::GetJob(int n, const GeantHPCWorker &worker) {
  GeantHPCJob job;
  job.type = JobType::Generator;
  job.events = std::min(n,eventsToDispatch);
  job.uid = GetNextId();
  eventsToDispatch -= job.events;
  return job;
}

bool GeantGeneratorJobPool::IsEmpty() {
  return eventsToDispatch == 0;
}

void GeantGeneratorJobPool::ReturnToPool(GeantHPCJob job) {
  eventsToDispatch += job.events;
}

void GeantGeneratorJobPool::SetEventAmount(int ev) {
  eventsToDispatch = ev;
}

void to_json(json& j, const GeantHepMCJob& p){
  j = json{{"fname",p.filename},{"offset",p.offset},{"n",p.amount}};
}

void from_json(const json& j, GeantHepMCJob& p){
  p.filename = j.at("fname").get<std::string>();
  p.amount = j.at("n").get<int>();
  p.offset = j.at("offset").get<int>();
}

GeantHPCJob GeantHPCJobPool::GetDublicateJob(GeantHPCJob &job) {
  GeantHPCJob newJob = job;
  newJob.uid = GetNextId();
  newJob.dublicateUIDs->insert(newJob.uid);
  return newJob;
}
}
}