
#include "TestFSData.h"

#include "HistFS.h"


namespace userfsapplication {

//
// TestFastSimDataPerPrimary
TestFastSimDataPerPrimary::TestFastSimDataPerPrimary() { Clear(); }

void TestFastSimDataPerPrimary::Clear() {
  fNumChargedSteps = fNumNeutralSteps = fChargedTrackL  = fNeutralTrackL = 0.;
  fNumGammas       = fNumElectrons    = fNumPositrons   = 0.;
  fNumPrimaryTrans = fNumPrimaryRefl  = fNumOneTrans    = fNumOneRefl = 0.;
  fEdepInTarget    = fELeakPrimary    = fELeakSecondary = 0.;
}

TestFastSimDataPerPrimary& TestFastSimDataPerPrimary::operator+=(const TestFastSimDataPerPrimary& other) {
  fNumChargedSteps += other.fNumChargedSteps;
  fNumNeutralSteps += other.fNumNeutralSteps;
  fChargedTrackL   += other.fChargedTrackL;
  fNeutralTrackL   += other.fNeutralTrackL;
  fNumGammas       += other.fNumGammas;
  fNumElectrons    += other.fNumElectrons;
  fNumPositrons    += other.fNumPositrons;
  fNumPrimaryTrans += other.fNumPrimaryTrans;
  fNumPrimaryRefl  += other.fNumPrimaryRefl;
  fNumOneTrans     += other.fNumOneTrans;
  fNumOneRefl      += other.fNumOneRefl;
  fEdepInTarget    += other.fEdepInTarget;
  fELeakPrimary    += other.fELeakPrimary;
  fELeakSecondary  += other.fELeakSecondary;
  return *this;
}




//
// TestFastSimData
TestFastSimData::TestFastSimData() { Clear(); }

void TestFastSimData::Clear() {
  fNumChargedSteps = fNumNeutralSteps = fNumChargedSteps2 = fNumNeutralSteps2 = 0.;
  fChargedTrackL   = fNeutralTrackL   = fChargedTrackL2   = fNeutralTrackL2   = 0.;
  fNumGammas       = fNumElectrons    = fNumPositrons     = 0.;
  fNumPrimaryTrans = fNumPrimaryRefl  = fNumOneTrans      = fNumOneRefl       = 0.;
  fEdepInTarget    = fEdepInTarget2   = 0.;
  fELeakPrimary    = fELeakSecondary  = fELeakPrimary2    = fELeakSecondary2  = 0.;
}

void TestFastSimData::AddDataPerPrimary(TestFastSimDataPerPrimary& data) {
  AddChargedSteps(data.GetChargedSteps());
  AddNeutralSteps(data.GetNeutralSteps());
  AddChargedTrackL(data.GetChargedTrackL());
  AddNeutralTrackL(data.GetNeutralTrackL());
  AddGammas   (data.GetGammas()   );
  AddElectrons(data.GetElectrons());
  AddPositrons(data.GetPositrons());
  AddPrimaryTransmitted(data.GetPrimaryTransmitted());
  AddPrimaryReflected(data.GetPrimaryReflected());
  AddOneTransmitted(data.GetOneTransmitted());
  AddOneReflected(data.GetOneReflected());
  AddEdepInTarget(data.GetEdepInTarget());
  AddELeakPrimary(data.GetELeakPrimary());
  AddELeakSecondary(data.GetELeakSecondary());
}




//
// TestFastSimDataPerEvent
TestFastSimDataPerEvent::TestFastSimDataPerEvent(int nprimperevent) : fNumPrimaryPerEvent(nprimperevent) {
  fPerPrimaryData.reserve(fNumPrimaryPerEvent);
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData.push_back(TestFastSimDataPerPrimary());
  }
}

void TestFastSimDataPerEvent::Clear() {
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData[i].Clear();
  }
}

TestFastSimDataPerEvent& TestFastSimDataPerEvent::operator+=(const TestFastSimDataPerEvent &other) {
  for (int i=0; i<fNumPrimaryPerEvent; ++i) {
    fPerPrimaryData[i] += other.fPerPrimaryData[i];
  }
  return *this;
}




//
// TestFastSimThreadDataEvents
TestFastSimThreadDataEvents::TestFastSimThreadDataEvents(int nevtbuffered, int nprimperevent) : fNumBufferedEvents(nevtbuffered) {
  fPerEventData.reserve(fNumBufferedEvents);
  for (int i=0; i<fNumBufferedEvents; ++i) {
    fPerEventData.push_back(TestFastSimDataPerEvent(nprimperevent));
  }
}

bool  TestFastSimThreadDataEvents::Merge(int evtslotindx, const TestFastSimThreadDataEvents& other) {
  fPerEventData[evtslotindx] += other.GetDataPerEvent(evtslotindx);
  return true;
}




//
// TestFastSimThreadDataRun
TestFastSimThreadDataRun::TestFastSimThreadDataRun() : fHisto1(nullptr) {}

TestFastSimThreadDataRun::~TestFastSimThreadDataRun() {
  if (fHisto1) {
    delete fHisto1;
  }
  fHisto1 = nullptr;
}

void TestFastSimThreadDataRun::CreateHisto1(int nbins, double min, double max) {
  if (fHisto1) {
    delete fHisto1;
  }
  fHisto1= new Hist(min, max, nbins);
}

bool TestFastSimThreadDataRun::Merge(int /*evtslotindx*/, const TestFastSimThreadDataRun& other) {
  (*fHisto1) += *(other.GetHisto1());
  return true;
}




} // namespace userfsapplication
