

#ifdef USE_ROOT
  #include "TH1F.h"
#else
  #include "Hist.h"
#endif

#include "ClicAppData.h"


namespace userapplication {

//
// ClicAppDataPerPrimary
ClicAppDataPerPrimary::ClicAppDataPerPrimary() { Clear(); }

void ClicAppDataPerPrimary::Clear() {

  for (int k=0;k<maxAbsorbers; k++){
  	fChargedTrackL[k]  = fNeutralTrackL[k] = fEdepInAbsorber[k] = 0.;
  }

  for (int k=0;k<numLayers; k++){
  	fChargedTrackLayerL[k]  = fNeutralTrackLayerL[k] = fEdepInLayer[k] = 0.;
  }


  fNumChargedSteps = fNumNeutralSteps = 0.;
  fNumGammas       = fNumElectrons    = fNumPositrons   = 0.;
  fELeakPrimary    = fELeakSecondary = 0.;
}

ClicAppDataPerPrimary& ClicAppDataPerPrimary::operator+=(const ClicAppDataPerPrimary& other) {
  for(int k=0;k<maxAbsorbers;k++) {
  	fChargedTrackL[k]   += other.fChargedTrackL[k];
  	fNeutralTrackL[k]   += other.fNeutralTrackL[k];
  	fEdepInAbsorber[k]    += other.fEdepInAbsorber[k];
  }
  for(int k=0;k<numLayers;k++) {
  	fChargedTrackLayerL[k]   += other.fChargedTrackLayerL[k];
  	fNeutralTrackLayerL[k]   += other.fNeutralTrackLayerL[k];
  	fEdepInLayer[k]    += other.fEdepInLayer[k];
  }
  fNumChargedSteps += other.fNumChargedSteps;
  fNumNeutralSteps += other.fNumNeutralSteps;
  fNumGammas       += other.fNumGammas;
  fNumElectrons    += other.fNumElectrons;
  fNumPositrons    += other.fNumPositrons;
  fELeakPrimary    += other.fELeakPrimary;
  fELeakSecondary  += other.fELeakSecondary;
  return *this;
}




//
// ClicAppData
ClicAppData::ClicAppData() { Clear(); }

void ClicAppData::Clear() {
  for (int k=1; k<maxAbsorbers;k++){
  	fChargedTrackL[k]   = fNeutralTrackL[k]   = fChargedTrackL2[k]   = fNeutralTrackL2[k]   = 0.;
  	fEdepInAbsorber[k]    = fEdepInAbsorber2[k]   = 0.;
  }

  fNumChargedSteps = fNumNeutralSteps = fNumChargedSteps2 = fNumNeutralSteps2 = 0.;
  fNumGammas       = fNumElectrons    = fNumPositrons     = 0.;
  fELeakPrimary    = fELeakSecondary  = fELeakPrimary2    = fELeakSecondary2  = 0.;
}

void ClicAppData::AddDataPerPrimary(ClicAppDataPerPrimary& data) {
  AddChargedSteps(data.GetChargedSteps());
  AddNeutralSteps(data.GetNeutralSteps());

  for (int k=0; k<maxAbsorbers; k++){
  	AddChargedTrackL(data.GetChargedTrackL(k),k);
  	AddNeutralTrackL(data.GetNeutralTrackL(k),k);
  	AddEdepInAbsorber(data.GetEdepInAbsorber(k),k);
  }

  for (int k=0; k<numLayers; k++){
  	AddChargedTrackLayerL(data.GetChargedTrackLayerL(k),k);
  	AddNeutralTrackLayerL(data.GetNeutralTrackLayerL(k),k);
  	AddEdepInLayer(data.GetEdepInLayer(k),k);
  }
  AddGammas   (data.GetGammas()   );
  AddElectrons(data.GetElectrons());
  AddPositrons(data.GetPositrons());


  AddELeakPrimary(data.GetELeakPrimary());
  AddELeakSecondary(data.GetELeakSecondary());
}

//
// ClicAppDataPerEvent
ClicAppDataPerEvent::ClicAppDataPerEvent(int nprimperevent) : fNumPrimaryPerEvent(nprimperevent) {
  fPerPrimaryData.reserve(fNumPrimaryPerEvent);
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData.push_back(ClicAppDataPerPrimary());
  }
}

void ClicAppDataPerEvent::Clear() {
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData[i].Clear();
  }
}

ClicAppDataPerEvent& ClicAppDataPerEvent::operator+=(const ClicAppDataPerEvent &other) {
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData[i] += other.fPerPrimaryData[i];
  }
  return *this;
}




//
// ClicAppDataEvents
ClicAppThreadDataEvents::ClicAppThreadDataEvents(int nevtbuffered, int nprimperevent) : fNumBufferedEvents(nevtbuffered) {
  fPerEventData.reserve(fNumBufferedEvents);
  for (int i=0; i<fNumBufferedEvents; ++i) {
    fPerEventData.push_back(ClicAppDataPerEvent(nprimperevent));
  }
}

bool ClicAppThreadDataEvents::Merge(int evtslotindx, const ClicAppThreadDataEvents& other) {
  fPerEventData[evtslotindx] += other.GetDataPerEvent(evtslotindx);
  return true;
}




//
// ClicAppThreadDataRun
ClicAppThreadDataRun::ClicAppThreadDataRun() : fHisto1(nullptr) {}

ClicAppThreadDataRun::~ClicAppThreadDataRun() {
  if (fHisto1) {
    delete fHisto1;
  }
  fHisto1 = nullptr;
}

void ClicAppThreadDataRun::CreateHisto1(int nbins, double min, double max) {
  if (fHisto1) {
    delete fHisto1;
  }
#ifdef USE_ROOT
  fHisto1= new TH1F("HistName", "Hist Title", nbins, min, max);
#else
  fHisto1= new Hist(min, max, nbins);
#endif
}

bool ClicAppThreadDataRun::Merge(int /*evtslotindx*/, const ClicAppThreadDataRun& other) {
#ifdef USE_ROOT
  fHisto1->Add(other.GetHisto1(),1);
#else
  (*fHisto1) += *(other.GetHisto1());
#endif
  return true;
}




} // namespace userapplication
