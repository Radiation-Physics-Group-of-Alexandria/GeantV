#include <HepMC/ReaderAscii.h>
#include <HepMC/ReaderRoot.h>
#include "HepMCGeneratorMultFiles.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

HepMCGeneratorMultFiles::HepMCGeneratorMultFiles() : filename("") , currentOffset(0) {

}
void HepMCGeneratorMultFiles::SetEventSource(const std::string &file, int offset) {
  if(filename != file || currentOffset > offset){ // we can't seek back
    filename = file;
    if(input_file) delete input_file;
    HepMCGenerator::LoadFile(filename);

    currentOffset = 0;
    desiredOffset = offset;
  }
}
GeantEventInfo HepMCGeneratorMultFiles::NextEvent() {
  if(currentOffset > desiredOffset)
    Geant::Fatal("HepMCGeneratorMultFiles::NextEvent", "Offset bigger that it should be.");

  HepMC::GenEvent dummyEvent(HepMC::Units::GEV, HepMC::Units::MM);
  while(currentOffset < desiredOffset){
    if (!(input_file->read_event(dummyEvent)))
      Geant::Fatal("HepMCGeneratorMultFiles::NextEvent", "No more events to read!");
    ++currentOffset;
  }

  return HepMCGenerator::NextEvent();
}

} // GEANT_IMPL_NAMESPACE
} // Geant