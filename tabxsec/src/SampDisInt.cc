//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: SampleInteractions.cc  J. Apostolakis  $
//
// -------------------------------------------------------------------
//      Geant4 Sampler file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      File name:     SampleInteractions
//
//      Author:        J. Apostolakis
//      Creation date: 20 June 2013
//
//      Adapted from:  test35.cc (HARP) by V.Ivanchenko
//
//      Modifications:
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4ElementVector.hh"

#include "G4VDiscreteProcess.hh"
#include "G4VContinuousDiscreteProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4ProcessManager.hh"
#include "G4FastStep.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleChangeForDecay.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4VRestProcess.hh"
#include "G4VEmProcess.hh"
#include "G4VRestDiscreteProcess.hh"

#include "G4ParticleTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Nucleus.hh"
#include "G4IonTable.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronStoppingProcess.hh"

#include "G4NucleiProperties.hh"

#include "G4UnitsTable.hh"
#include "G4StateManager.hh"

#include "G4NistManager.hh"

#include "G4ForceCondition.hh"
#include "G4TouchableHistory.hh"

#include "G4VCrossSectionDataSet.hh"

#include "G4TransportationManager.hh"
#include "G4Navigator.hh"

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "G4Timer.hh"

#include "SampleDisInt.hh"

#include <TFinState.h>

using namespace std;

int SampDisInt(G4Material *material, G4ThreeVector *pos, G4DynamicParticle *dpart, G4double de, G4VProcess *proc,
               G4int nevt, G4int verbose, TFinState &fs) {

  G4VDiscreteProcess *discProc = dynamic_cast<G4VDiscreteProcess *>(proc);
  G4VContinuousDiscreteProcess *contdProc = dynamic_cast<G4VContinuousDiscreteProcess *>(proc);
  G4VRestDiscreteProcess *discRestProc = dynamic_cast<G4VRestDiscreteProcess *>(proc);

  const G4String pname = proc->GetProcessName();

  if (!discProc && !contdProc && !discRestProc) {
    G4cout << " Process " << pname << " for particle " << dpart->GetParticleDefinition()->GetParticleName()
           << " has no discrete part. " << G4endl;
    return -1;
  }

  if (G4String("Decay") == pname) {
    G4DecayTable *dt = dpart->GetParticleDefinition()->GetDecayTable();
    if (!dt) {
      G4cout << "Cannot decay " << dpart->GetParticleDefinition()->GetParticleName() << " without a decay table"
             << G4endl;
      return -1;
    }
  }

  // ------- Define target A
  const G4Element *elm = material->GetElement(0);
  G4int A = (G4int)(elm->GetN() + 0.5);
  G4int Z = (G4int)(elm->GetZ() + 0.5);
  G4double amass = GetNuclearMass(Z, A, verbose);

  Finstat_t *fstat = new Finstat_t[nevt];
  memset(fstat, 0, nevt * sizeof(Finstat_t));

  if (verbose > 1) {
    G4cout << G4endl << "Process   : " << pname << G4endl
           << "Particle  : " << dpart->GetParticleDefinition()->GetParticleName() << "; p " << dpart->Get4Momentum()
           << G4endl << "Material  : " << material->GetName() << "  Z " << Z << "  A " << A << " mass " << amass
           << G4endl << " Position " << *pos / mm << "mm" << G4endl;
  }

  //----------------------------- Start event loop -------------------------------
  G4int iter = 0;
  G4int nresample = 0;
  const G4int maxresample = 10 * nevt;
  G4double eini = dpart->GetKineticEnergy();
  while (iter < nevt) {
    if (de != 0) {
      G4double ep = eini + de * G4UniformRand();
      dpart->SetKineticEnergy(ep);
    }
    if (verbose > 1) {
      G4cout << "### " << iter << "-th event: " << dpart->GetParticleDefinition()->GetParticleName() << " " << pname
             << " @" << dpart->Get4Momentum() << " on " << material->GetName() << " (" << Z << "," << A << "," << amass
             << ")";
      // ------- Printout
    }

    if (++nresample > maxresample)
      break;
    SampleOne(material, pos, dpart, proc, verbose, fstat[iter]);
    if (pname != G4String("Decay") || fstat[iter].npart) {
      ++iter;
    } else {
      G4cout << "No particles produced for decay!!!!" << G4endl;
    }
  }

  //----------------------------- End event loop -------------------------------

  //----------------------------- Build final structure ------------------------
  G4int ntotp = 0;
  for (G4int i = 0; i < nevt; ++i)
    ntotp += fstat[i].npart;

  //  printf("Total particles %d\n",ntotp);

  if (ntotp) {
    G4int *npart = new G4int[nevt];
    G4float *weight = new G4float[nevt];
    G4float *kerma = new G4float[nevt];
    G4float *en = new G4float[nevt];
    char *surv = new char[nevt];

    G4int *tpid = new G4int[ntotp];
    G4float *tmom = new G4float[3 * ntotp];

    G4int ip = 0;
    for (G4int i = 0; i < nevt; ++i) {
      npart[i] = fstat[i].npart;
      weight[i] = fstat[i].weight;
      kerma[i] = fstat[i].kerma / GeV;
      en[i] = fstat[i].en / GeV;
      surv[i] = fstat[i].survived;
      for (G4int j = 0; j < npart[i]; ++j) {
        tpid[ip + j] = fstat[i].pid[j];
        tmom[3 * (ip + j)] = fstat[i].mom[3 * j] / GeV;
        tmom[3 * (ip + j) + 1] = fstat[i].mom[3 * j + 1] / GeV;
        tmom[3 * (ip + j) + 2] = fstat[i].mom[3 * j + 2] / GeV;
      }
      delete[] fstat[i].pid;
      delete[] fstat[i].mom;
      ip += npart[i];
    }

    fs.SetFinState(nevt, npart, weight, kerma, en, surv, tpid, tmom);
    //    fs.Print();

    delete[] npart;
    delete[] weight;
    delete[] kerma;
    delete[] en;
    delete[] surv;

    delete[] tmom;
    delete[] tpid;
  }

  delete[] fstat;
  return true;
}

void TimingInfo(G4float, G4int, G4int, G4int, G4float, G4int);

G4int SampleOne(G4Material *material, G4ThreeVector *pos, G4DynamicParticle *dpart, G4VProcess *proc, G4int verbose,
                Finstat_t &fs)

{

  // G4VDiscreteProcess *discProc = dynamic_cast<G4VDiscreteProcess *>(proc);
  G4VContinuousDiscreteProcess *contdProc = dynamic_cast<G4VContinuousDiscreteProcess *>(proc);
  G4HadronicProcess *hadp = dynamic_cast<G4HadronicProcess *>(proc);
  G4VEmProcess *vemp = dynamic_cast<G4VEmProcess *>(proc);
  G4VRestProcess *restProc = dynamic_cast<G4VRestProcess *>(proc);
  G4VRestDiscreteProcess *restDiscProc = dynamic_cast<G4VRestDiscreteProcess *>(proc); // e.g. decay
  G4HadronStoppingProcess *hStoppingProc = dynamic_cast<G4HadronStoppingProcess *>(proc);
  const G4String pname = proc->GetProcessName();

  // Timing stuff
  clock_t begin, end;
  //

  // -------- We chose a short step because we do not want anything
  // -------- to happen during the step

  const G4double theStep = 1 * mm;

  // -------- Projectile
  const G4ThreeVector &aPosition = *pos;
  G4ThreeVector aDirection = dpart->GetMomentumDirection();
  const G4DynamicParticle *dpsave = new G4DynamicParticle(*dpart);

  G4double mass = dpart->GetParticleDefinition()->GetPDGMass();
  G4String name = dpart->GetParticleDefinition()->GetParticleName();

  G4double e0 = dpart->GetKineticEnergy();
  G4double charge = dpart->GetParticleDefinition()->GetPDGCharge(); // Baseline charge
  // During tracking use dpart->GetCharge() instead - as the charge can be dynamic (e.g. for ions)

  // G4VCrossSectionDataSet* GetCrossSectionDS(const G4ParticleDefinition*, G4Material *);
  // G4VCrossSectionDataSet* cs = GetCrossSectionDS(part, material);

  // -------- Track
  G4Track *gTrack = new G4Track(new G4DynamicParticle(*dpart), 0, G4ThreeVector(aPosition));
  G4TouchableHandle touchable = new G4TouchableHistory();
  G4TransportationManager::GetTransportationManager()
      ->GetNavigatorForTracking()
      ->LocateGlobalPointAndUpdateTouchableHandle(aPosition, aDirection, touchable, false);

  gTrack->SetTouchableHandle(touchable);

  // -------- Step

  G4Step *step;
  step = new G4Step();

  // Intertwines the step and track - does not Set their values
  step->SetTrack(gTrack);
  gTrack->SetStep(step);

  G4StepPoint *aPoint, *bPoint;
  aPoint = new G4StepPoint();
  aPoint->SetPosition(aPosition);
  aPoint->SetMaterial(material);
  G4double safety = 10000. * cm;
  aPoint->SetSafety(safety);
  aPoint->SetTouchableHandle(touchable);

  aPoint->SetKineticEnergy(e0);
  aPoint->SetMass(mass); //  To be sure
  aPoint->SetCharge(charge);

  step->SetPreStepPoint(aPoint);

  const G4MaterialCutsCouple *mcc = FindMaterialCutsCouple(material);

  aPoint->SetMaterialCutsCouple(mcc);
  aPoint->SetMaterial(material);

  bPoint = new G4StepPoint(*aPoint);
  G4ThreeVector bPosition = aPosition + aDirection * theStep;
  bPoint->SetPosition(bPosition);
  step->SetPostStepPoint(bPoint);
  step->SetStepLength(theStep);

  if (!G4StateManager::GetStateManager()->SetNewState(G4State_Idle))
    cout << "G4StateManager PROBLEM: Not able to set it to Idle state! " << G4endl;

  // G4double e;
  G4VParticleChange *aChange = 0;

  G4int isurv = 0;
  fs.survived = FALSE;
  G4ThreeVector psurv(0, 0, 0);
  G4int spid = 0;

  fs.kerma = 0;
  G4bool needendl = FALSE;
  // G4int modu = 10000;

  // gTrack->SetKineticEnergy(e0);
  if (std::fabs(gTrack->GetKineticEnergy() - e0) > 1.0e-8 * e0) {
    G4cerr << " ERROR> : Track Kinetic energy is not set correctly " << G4endl;
    G4cerr << "        :  Track Kin E = " << gTrack->GetKineticEnergy() << G4endl;
    G4cerr << "        :  Expected    = " << e0 << G4endl;
    exit(1);
    // (gTrack->GetKineticEnergy - e0)
  }
  gTrack->SetStep(step);
  proc->StartTracking(gTrack);

  G4LorentzVector labv;
  labv = G4LorentzVector(gTrack->GetMomentum(), gTrack->GetTotalEnergy());
  if (labv != dpsave->Get4Momentum()) {
    if (needendl) {
      G4cout << G4endl;
      needendl = FALSE;
    }
    G4cout << "The momentum of the projectile has changed from " << dpsave->Get4Momentum() << " for " << name << " to "
           << labv << " for " << gTrack->GetParticleDefinition()->GetParticleName() << G4endl;
    G4cout << gTrack->GetMomentum() << gTrack->GetTotalEnergy() << gTrack->GetParticleDefinition()->GetPDGMass()
           << G4endl;
    exit(1);
  }

  // ----------------------------------- ** Cause Discrete Interaction ** ----------------------------

  begin = clock();

  if (e0 == 0.0 && (restProc || restDiscProc)) {
    // G4ForceCondition fCondition;
    // track status must be set properly
    gTrack->SetTrackStatus(fStopButAlive);
    // In SteppingManager::DefinePhysicalStepLength()
    // G4double physIntLength = proc->AtRestGetPhysicalInteractionLength(*gTrack, &fCondition);
    // -- Make it happen in time
    aChange = proc->AtRestDoIt(*gTrack, *step);
    aChange->UpdateStepForAtRest(step);

  } else {

    G4double previousStepSize = 1.0;
    G4ForceCondition fCondition;
    //G4double physIntLength;

    // In SteppingManager::DefinePhysicalStepLength()
    if (proc->isPostStepDoItIsEnabled())
      /*physIntLength = */proc->PostStepGetPhysicalInteractionLength(*gTrack, previousStepSize, &fCondition);

    // Ignore the proposed value - will force the interaction
    //  fCondition= Forced;

    if (contdProc) {
      // -- if continuous discrete process, make it happen along the step
      G4GPILSelection fGPILSelection;
      G4double safetyPrx = 2 * theStep; // Not limiting
      /*physIntLength = */contdProc->AlongStepGetPhysicalInteractionLength(*gTrack, previousStepSize, theStep, safetyPrx,
                                                                       &fGPILSelection);

      // safetyPrx is output only for transportation - currently

      aChange = contdProc->AlongStepDoIt(*gTrack, *step);
      // Update the PostStepPoint of Step according to ParticleChange
      aChange->UpdateStepForAlongStep(step);

      // ignore any energy loss along the step: so we don't update gTrack and
      // reset Edepo and clear particleChange before the discrete part. Along
      // step energy loss will be computed in the MC from dE/dx
      step->ResetTotalEnergyDeposit();
      aChange->Clear();

      // Check the energy loss during the Along Step
      if (std::fabs(gTrack->GetKineticEnergy() - e0) > 0.5 * e0) {
        G4cerr << " Warning : Track lost majority of its Kinetic energy during AlongStep " << G4endl;
        G4cerr << "      Current Kin E (T) = " << gTrack->GetKineticEnergy() << G4endl;
        G4cerr << "      Initial Kin E (T0)= " << e0 << G4endl;
        G4cerr << "      Ratio (T/T0)  = " << gTrack->GetKineticEnergy() / e0 << G4endl;
      }
    }

    // need to do this because of the inheritance problem in G4HadronicProcess
    if (proc->isPostStepDoItIsEnabled()) {
      // -- Make it happen at the end of the step
      aChange = proc->PostStepDoIt(*gTrack, *step);
      // std::cout << "SDI> PostStepDoIt: Nsec= " << aChange->GetNumberOfSecondaries() << std::endl;
      aChange->UpdateStepForPostStep(step);
    } else if (e0 == 0.0 && proc->isAtRestDoItIsEnabled()) {
      gTrack->SetTrackStatus(fStopButAlive);
      /*physIntLength = */proc->AtRestGetPhysicalInteractionLength(*gTrack, &fCondition);
      aChange = proc->AtRestDoIt(*gTrack, *step);
      aChange->UpdateStepForAtRest(step);
    }
  }
  step->UpdateTrack();

  end = clock();
  G4float cputime = (G4double)(end - begin) / CLOCKS_PER_SEC;

  G4int n = aChange->GetNumberOfSecondaries();

  fs.kerma += step->GetTotalEnergyDeposit();

  // ------- Define target A
  const G4Element *elm = material->GetElement(0);
  G4int A = (G4int)(elm->GetN() + 0.5);
  G4int Z = (G4int)(elm->GetZ() + 0.5);
  G4double amass = GetNuclearMass(Z, A, verbose);
  G4int isoA = A;
  G4int isoZ = Z;
  G4double isoM = amass;
  if (hadp) {
    const G4Isotope *iso = hadp->GetTargetIsotope();
    if (!iso) {
      if (verbose)
        G4cout << "Process " << proc->GetProcessName() << " did not select isotope!" << G4endl;
    } else {
      //--this will results in a mass in which the mass of the electron shell is
      // included as well. But we want the mass of the target nucleus so replace
      // this to get it correctly. This will resolve some energy cons. violations
      // obtained earlier.
      // isoM = iso->GetA()/Avogadro*c_squared;
      isoZ = iso->GetZ();
      isoA = iso->GetN();
      isoM = G4NucleiProperties::GetNuclearMass(isoA, isoZ);
      if (verbose)
        if ((A != isoA) || (Z != isoZ))
          printf("Target changed. A: %d -> %d, Z: %d -> %d, M: %f -> %f (%f)\n", A, isoA, Z, isoZ, amass, isoM,
                 GetNuclearMass(isoZ, isoA, verbose));
    }
  }

  if (vemp) {
    // --- "by hand" determination of the target isotope, it is not saved in G4
    const G4Element *ele = vemp->GetCurrentElement();
    if (!ele) {
      if (verbose > 1)
        G4cout << "Process " << vemp->GetProcessName() << " did not select isotope" << G4endl;
    } else {
      isoZ = ele->GetZ();
      if (isoZ == Z) {
        // the target nucleus is what we think it should be
        for (G4int i = 0; i < n; ++i) {
          G4Ions *prion = dynamic_cast<G4Ions *>(aChange->GetSecondary(i)->GetDefinition());
          if (prion) {
            if (prion->GetPDGCharge() == Z) {
              isoA = prion->GetBaryonNumber();
              isoM = GetNuclearMass(isoZ, isoA, verbose);
              break;
            }
          }
        }
      }
      if (verbose > 1)
        G4cout << "Interaction " << vemp->GetProcessName() << " happened on Z " << ele->GetZ() << " A " << ele->GetN()
               << G4endl;
    }
  }

  // *** Warning *** this is experimental *** We use the Target Isotope for energy and B balance

  A = isoA;
  amass = isoM;

  // ----------------------------------- Correct for checking energy balance --------------

  G4LorentzVector pcons(labv);

  // -- add the mass and baryon number of the target in case it is a hadron inelastic
  G4int bnum = dpart->GetParticleDefinition()->GetBaryonNumber();
  if (dynamic_cast<G4HadronInelasticProcess *>(proc) || dynamic_cast<G4VRestProcess *>(proc) ||
      dynamic_cast<G4HadronCaptureProcess *>(proc)) {
    bnum += A;
    pcons[3] += amass;
  }

  // -- add the mass of the electron in case it is an annihilation
  if (name == G4String("e+") && pname == G4String("annihil"))
    pcons[3] += mass;

  //
  // We save at this point the value of the "input energy" to check the energy balance
  // This is highly debatable, but it is a good compromise

  G4LorentzVector porig = pcons;

  // -- See wheter the original particle is still alive
  if ((aChange->GetTrackStatus() == fAlive) || (aChange->GetTrackStatus() == fStopButAlive)) {
    // -- Original particle is alive -- let's get it
    if (verbose > 1)
      printf("after %s on %s %s is still alive\n", (const char *)pname, (const char *)material->GetName(),
             (const char *)name);

    // printf("KERMA 1 := %f kinEnergy 1:= %f\n", fs.kerma, gTrack->GetKineticEnergy());
    // !step->GetTrack()->GetMomentum() is a null-vector if the status is fStop...
    // because kinetic energy of the particle is set to zero. Momentum conservation
    // cannot be checked in these cases and will be skipped.!
    // we remove the life particle from the input vector to test momentum conservation
    pcons -= G4LorentzVector(step->GetTrack()->GetMomentum(), step->GetTrack()->GetTotalEnergy());
    bnum -= dpart->GetParticleDefinition()->GetBaryonNumber();

    fs.survived = TRUE;
    isurv = 1;
    psurv = step->GetTrack()->GetMomentum();
    spid = TPartIndex::I()->PartIndex(dpart->GetParticleDefinition()->GetPDGEncoding());
    if (G4String("hadElastic") == pname) {
      // ----------------------------------- Trying to find out the angle in elastic ---------------------------
      //--this part is meaningless if the status is not fAlive because
      // gTrack->GetMomentum() is a null-vector and gTrack->GetKineticEnergy() is
      // zero if the status is fStop... . So do it only
      if (aChange->GetTrackStatus() == fAlive) {
        G4double cost = labv.vect().cosTheta(gTrack->GetMomentum());
        // G4double ken = gTrack->GetKineticEnergy();
        G4double mas = gTrack->GetParticleDefinition()->GetPDGMass();
        // G4double pmo = gTrack->GetMomentum().mag();
        if (needendl) {
          G4cout << G4endl;
          needendl = FALSE;
        }
        if (verbose)
          G4cout << "Elastic scattering by cosTheta " << cost << " (" << 180 * std::acos(cost) / std::acos(-1)
                 << " deg)"
                 << " Output p " << G4LorentzVector(gTrack->GetMomentum(), gTrack->GetTotalEnergy()) << " mass " << mas
                 << G4endl;
      }

      // we add the baryon number of the target but only if there are  generated particles
      // oterwise the target will recoil and this will alter the baryon number conservation
      if (n) {
        bnum += A;
        pcons[3] += amass;
        porig[3] += amass;
      }
    }
  }

  if (G4String("CoulombScat") == pname) {
    // If the target nucleus recoils we have to add the mass to the initial balance
    if (verbose)
      G4cout << "CoulombScatt status " << tStatus[aChange->GetTrackStatus()] << " p "
             << G4LorentzVector(step->GetTrack()->GetMomentum(), step->GetTrack()->GetTotalEnergy()) << G4endl;
    if (n) {
      pcons[3] += amass;
      porig[3] += amass;
      bnum += A;
    }
  }

  TimingInfo(cputime, Z, dpart->GetParticleDefinition()->GetPDGEncoding(),
             proc->GetProcessType() * 1000 + proc->GetProcessSubType(), dpart->GetKineticEnergy(), n + isurv);

  // This is ok. The particle energy drops below 'theLowestEnergy' (set to 1keV
  // in G4HardonElasticProcess) due to recoil effect and the particle doesn't
  // have any AtRest processes. So its kinetic energy and status are set to zero
  // and fStopAndKill.
  /*
  if(G4String("hadElastic") == pname && !fs.survived) {
    G4cout << "Elastic but the particle did not survive " << tStatus[aChange->GetTrackStatus()] << G4endl;

  }
  */

  if (n)
    if ((G4String("hIoni") == pname) || (G4String("ionIoni") == pname) || (G4String("eIoni") == pname) ||
        (G4String("phot") == pname) || (G4String("compt") == pname)) {
      // the electron is NOT produced in the collision ... ;-)
      pcons[3] += G4ParticleTable::GetParticleTable()->FindParticle("e-")->GetPDGMass();
      porig[3] += G4ParticleTable::GetParticleTable()->FindParticle("e-")->GetPDGMass();
    }

  G4int iextra = 0;
  if (proc->GetProcessType() == 4 && proc->GetProcessSubType() == 131) {
    // Troubles with capture: the resulting isotope is not generated by G4, so there is imbalance in
    // momentum and baryon number. We try to correct here, but it is going to be messy...
    // Make sure that we do it for neutrons only at the beginning
    if (G4String("neutron") == name) {
      // first make sure that G4 did not put a ion in the output particles
      G4bool lprion = FALSE;
      for (G4int i = 0; i < n; ++i) {
        G4Ions *prion = dynamic_cast<G4Ions *>(aChange->GetSecondary(i)->GetDefinition());
        if (prion) {
          lprion = TRUE;
          if (prion->GetPDGCharge() != Z || prion->GetBaryonNumber() != A + 1) {
            // G4 put an output ion but it does not seem right. Put an error message and continue
            G4cout << name << " " << pname << " on " << material->GetName() << " (" << Z << "," << A << "," << amass
                   << ")"
                   << " produces Z " << prion->GetPDGCharge() << " A " << prion->GetBaryonNumber() << G4endl;
          }
        }
      }
      if (!lprion) {
        // No ion was produced by G4, so we can add our own
        iextra = 1;
      }
      // G4cout << "Capture is " << proc->GetProcessName() << G4endl;
    }
  }

  G4int ntot = isurv + n + iextra;
  fs.npart = 0;
  fs.weight = 1; // for the moment all events have the same prob
  fs.en = gTrack->GetKineticEnergy();
  // printf("KERMA 2 := %f kinEnergy 2:= %f\n", fs.kerma, fs.en);

  if (ntot) {
    fs.mom = new G4float[3 * (ntot)];
    fs.pid = new G4int[ntot];
    if (isurv) {
      fs.mom[0] = psurv[0];
      fs.mom[1] = psurv[1];
      fs.mom[2] = psurv[2];
      fs.pid[0] = spid;
      fs.npart = 1;
    }
  } else {
    fs.mom = 0;
    fs.pid = 0;
  }

  G4int nbar = 0;
  G4int n_pr = 0;
  G4int n_nt = 0;
  G4int n_pi = 0;
  G4int n_ga = 0;
  G4int n_el = 0;
  G4int n_po = 0;
  const G4DynamicParticle *sec = 0;
  G4ParticleDefinition *pd;

  //  Examine the secondaries
  G4DynamicParticle *secs = 0;

  // ----------------------------------- Generated secondaries ----------------------------
  // post-interaction primary
  if (verbose > 1) {
    G4cout << "POST INTERACTION INFO: " << G4endl;
    G4cout << "E-depo: gTrack->GetStep()->GetTotalEnergyDeposit() = " << gTrack->GetStep()->GetTotalEnergyDeposit()
           << G4endl;
    if (gTrack->GetTrackStatus() == fAlive || gTrack->GetTrackStatus() == fStopButAlive) {
      sec = gTrack->GetDynamicParticle();
      pd = sec->GetDefinition();
      G4cout << "Post-inter. primary " << pd->GetParticleName() << " (" << pd->GetPDGEncoding() << ") "
             << " Z= " << pd->GetAtomicNumber() << " B= " << pd->GetBaryonNumber() << " p= " << sec->Get4Momentum()
             << " Ekin =  " << sec->GetKineticEnergy() << G4endl;
    }
  }

  if (n)
    secs = new G4DynamicParticle[n];

  for (G4int i = 0; i < n; ++i) {

    sec = aChange->GetSecondary(i)->GetDynamicParticle();
    secs[i] = *sec;
    pd = sec->GetDefinition();
    G4int enc = pd->GetPDGEncoding();
    if (enc == -2212)
      ++nbar;
    if (enc == 2212)
      ++n_pr;
    if (enc == 2112)
      ++n_nt;
    if (enc == 22)
      ++n_ga;
    if (enc == 11)
      ++n_el;
    if (enc == -11)
      ++n_po;
    if (std::fabs(enc) == 211 || enc == 210) {
      ++n_pi;
    }

    if (verbose > 1) {
      if (needendl) {
        G4cout << G4endl;
        needendl = FALSE;
      }
      G4cout << " Sec[" << i << "]=" << pd->GetParticleName() << " (" << pd->GetPDGEncoding() << ") "
             << " Z= " << pd->GetAtomicNumber() << " B= " << pd->GetBaryonNumber() << " p= " << sec->Get4Momentum()
             << " Ekin =  " << sec->GetKineticEnergy() << G4endl;
    }
  }

  if (verbose > 1) {
    G4cout << ": " << n << " sec (" << n_pr << " protons, " << n_nt << " neutrons, " << n_pi << " pi), " << n_el
           << " e-, " << n_po << " e+, " << n_ga << " g" << G4endl;
  }

  if (proc->GetProcessType() == 4 && proc->GetProcessSubType() == 131 && iextra) {
    // This is a capture and we have to add our own output isotope balancing momentum - energy
    G4LorentzVector ptemp = pcons;
    for (G4int i = 0; i < n; ++i) {
      ptemp -= secs[i].Get4Momentum();
    }
    G4double isomass = GetNuclearMass(Z, A + 1, verbose);
    // G4double mmass = ptemp.m();
    ptemp[3] = std::sqrt(ptemp.vect().mag2() + isomass * isomass);
    pcons -= ptemp;
    bnum -= A + 1;
    if (verbose)
      G4cout << "nCapture: Mass of missing momentum " << ptemp.m() << " should be " << isomass << G4endl;
    // Adding output ion
    G4int ip = fs.npart;
    fs.pid[ip] = 1000000000 + 10000 * Z + 10 * (A + 1);
    fs.mom[3 * ip] = ptemp[0];
    fs.mom[3 * ip + 1] = ptemp[1];
    fs.mom[3 * ip + 2] = ptemp[2];
    ++fs.npart;
  }

  if (n) {                                        // We check only if we have secondaries
    if ((G4String("conv") != pname)               // But not if we have conversion, momentum is not preserved
        && (G4String("PositronNuclear") != pname) // Here we have a problem in G4
        && (G4String("ElectroNuclear") != pname)  // Another problem in G4
        &&
        !((G4String("hadElastic") == pname) &&
          (aChange->GetTrackStatus() != fAlive)) // post interaction momentum is lost
        &&
        !hStoppingProc // would be complicated, so skipp now
        ) {
      G4double perr;
      G4int berr;
      G4LorentzVector ptest;

      const G4double prec = 1e-2;
      checkBalance(porig, pcons, bnum, secs, n, ptest, perr, berr);

      if (perr > prec || berr) {
        // --- Try to fix few simple things
        if (ptest[3] > -939 && ptest[3] < -936 && berr == -1) {
          // G4cout << "Energy balance wrong by one baryon mass, let's try to fix it changing a deuteron in proton" <<
          // G4endl;
          for (G4int i = 0; i < n; ++i) {
            if (G4String("deuteron") == secs[i].GetDefinition()->GetParticleName()) {
              secs[i] = G4DynamicParticle(G4Proton::Proton(), secs[i].GetMomentum());
              break;
            }
          }
          checkBalance(porig, pcons, bnum, secs, n, ptest, perr, berr);
          if (ptest[3] > -939 && ptest[3] < -936 && berr == -1) {
            // -- this still did not fix it... however sometimes there is a neutron that we can change into a gamma
            // G4cout << "Energy balance wrong by one baryon mass, let's try to fix it changing a neutron in gamma" <<
            // G4endl;
            for (G4int i = 0; i < n; ++i) {
              if (G4String("neutron") == secs[i].GetDefinition()->GetParticleName()) {
                secs[i] = G4DynamicParticle(G4Gamma::Gamma(), secs[i].GetMomentum());
                break;
              }
            }
          }
        }
        checkBalance(porig, pcons, bnum, secs, n, ptest, perr, berr);
        if (perr > prec || berr) {
          pd = gTrack->GetDefinition();
          G4cout << setfill('-') << setw(120) << "-" << setfill(' ') << setw(0) << G4endl << "Dubious E/p/B balance "
                 << setiosflags(ios::scientific) << setprecision(2) << perr << " / dB " << berr
                 << " delta p=" << setiosflags(ios::fixed) << setprecision(6) << ptest << " for " << pname << " of "
                 << name << " @ " << dpart->Get4Momentum() << " on " << material->GetName() << " (" << Z << "," << A
                 << "," << amass << ")" << G4endl
                 //        << "  Out   :"
                 //        << setiosflags(ios::left) << setw(10) << pd->GetParticleName()
                 //        << " (" << setiosflags(ios::right) << setw(6) << pd->GetPDGEncoding() << ") " << setw(0)
                 //        << " Z:" << pd->GetAtomicNumber() << " B:" << pd->GetBaryonNumber()
                 //        << " p:" << gTrack->GetDynamicParticle()->Get4Momentum() << ": " <<
                 //        tStatus[aChange->GetTrackStatus()] << G4endl

                 << "  Out   :" << setiosflags(ios::left) << setw(10)
                 << step->GetTrack()->GetParticleDefinition()->GetParticleName() << " (" << setiosflags(ios::right)
                 << setw(6) << step->GetTrack()->GetParticleDefinition()->GetPDGEncoding() << ") " << setw(0)
                 << " Z:" << step->GetTrack()->GetParticleDefinition()->GetAtomicNumber()
                 << " B:" << step->GetTrack()->GetParticleDefinition()->GetBaryonNumber()
                 << " p:" << G4LorentzVector(step->GetTrack()->GetMomentum(), step->GetTrack()->GetTotalEnergy())
                 << ": " << tStatus[aChange->GetTrackStatus()] << G4endl

                 << "Particles generated:" << G4endl;
          for (G4int i = 0; i < n; ++i) {

            sec = &secs[i];
            pd = sec->GetDefinition();
            G4cout << "  #" << setw(3) << i << setw(0) << ": " << setiosflags(ios::left) << setw(10)
                   << pd->GetParticleName() << " (" << setiosflags(ios::right) << setw(6) << pd->GetPDGEncoding()
                   << ") " << setw(0) << " Z:" << pd->GetAtomicNumber() << " B:" << pd->GetBaryonNumber()
                   << " p:" << sec->Get4Momentum() << G4endl;
          }
          G4cout << setfill('-') << setw(120) << "-" << setfill(' ') << setw(0) << G4endl;
        }
      }
    }
    for (G4int i = 0; i < n; ++i) {
      G4ParticleDefinition *pde = secs[i].GetDefinition();
      G4int enc = pde->GetPDGEncoding();
      G4int gVpid = 0;
      if (enc > 1000000000 && enc < 2000000000)
        gVpid = enc;
      else
        gVpid = TPartIndex::I()->PartIndex(enc);
      if (gVpid < 0) {
        // If there is Ion we put it in the output stack. However we have to handle it afterwards...
        fs.kerma += secs[i].GetKineticEnergy();
      } else {
        G4int ip = fs.npart;
        fs.pid[ip] = gVpid;
        fs.mom[3 * ip] = secs[i].Get4Momentum()[0];
        fs.mom[3 * ip + 1] = secs[i].Get4Momentum()[1];
        fs.mom[3 * ip + 2] = secs[i].Get4Momentum()[2];
        ++fs.npart;
      }
    }
  }

  for (G4int i = 0; i < n; ++i)
    delete aChange->GetSecondary(i);
  aChange->Clear();

  delete[] secs;

  // A step owns its Step points
  //   - must ensure that they are valid or null (and not pointint to same StepPt)
  delete step;

  // A track owns its dynamic particle - that will be deleted by it
  delete gTrack;

  return 1;
}

G4double GetNuclearMass(G4int Z, G4int N, G4int verbose) {
  G4double mass = 0.0;
  G4Nucleus targetNucleus;

  if (verbose > 1)
    G4cout << "Nucleus with N= " << N << "  Z= " << Z << G4endl;
  targetNucleus.SetParameters((G4double)N, (G4double)Z);
  mass = targetNucleus.AtomicMass((G4double)N, (G4double)Z);
  if (verbose > 1)
    G4cout << "Mass from targetNucleus(MeV)= " << mass / MeV << G4endl;
  mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, N);
  if (verbose > 1)
    G4cout << "Mass from IonTable(MeV)=      " << mass / MeV << G4endl;
  return mass;
}

#include "G4ProductionCutsTable.hh"
// G4MaterialCutsCouple

const G4MaterialCutsCouple *FindMaterialCutsCouple(G4Material *mat) {
  const G4MaterialCutsCouple *couple = 0;

  G4ProductionCutsTable *theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();

  size_t numOfCouples = theCoupleTable->GetTableSize();

  for (size_t i = 0; i < numOfCouples; i++) {
    couple = theCoupleTable->GetMaterialCutsCouple(i);
    if (couple->GetMaterial() == mat)
      break;
  }
  if (couple == 0) {
    G4cerr << "Fatal ERROR> Not found couple for material " << mat->GetName() << G4endl;
    exit(1);
  }
  return couple;
}

void checkBalance(const G4LorentzVector &porig, const G4LorentzVector &pmom, G4int bnum, G4DynamicParticle *secs,
                  G4int n, G4LorentzVector &ptest, G4double &perr, G4int &berr) {
  ptest = pmom;
  berr = bnum;
  for (G4int i = 0; i < n; ++i) {
    ptest -= secs[i].Get4Momentum();
    berr -= secs[i].GetDefinition()->GetBaryonNumber();
  }

  G4double err = 0;
  G4double ptot = porig.vect().mag();
  perr = std::abs(ptest[3]);
  if (pmom[3])
    perr /= porig[3];
  if (ptot)
    for (G4int i = 0; i < 3; ++i) {
      err = std::abs(ptest[i]) / ptot;
      if (err > perr)
        perr = err;
    }
}

G4bool rescaleEnergy(const G4LorentzVector &porig, G4DynamicParticle *secs, G4int n, G4double eleft, G4double etot) {
  //  G4double esum = etot;
  G4LorentzVector psum(0, 0, 0, 0);
  for (G4int i = 0; i < n; ++i)
    psum += secs[i].Get4Momentum();
  if (std::abs(eleft - (etot - psum[3])) > 1e-6) {
    G4cout << "You screwed it up" << G4endl;
    exit(1);
  }
  //  G4DynamicParticle *bsecs = new G4DynamicParticle[n];
  // The cms momentum
  G4ThreeVector cms(psum[0] / psum[3], psum[1] / psum[3], psum[2] / psum[3]);
  G4ThreeVector cmsori(porig[0] / porig[3], porig[1] / porig[3], porig[2] / porig[3]);

  // -- see what is the total energy if there is just no motion
  G4LorentzVector pcms(0, 0, 0, 0);
  for (G4int i = 0; i < n; ++i)
    pcms[3] += secs[i].GetDefinition()->GetPDGMass();
  G4LorentzVector plab = pcms.boost(cms);
  G4cout << "Energy of all particle at rest " << plab << " Energy available " << porig[3] << " Delta "
         << plab[3] - porig[3] << G4endl;

  if (plab[3] < porig[3]) {
    // The energy of the bare masses is less than the energy of the reaction, we can therefore rescale
    // in both directions
    G4cout << "Here rescaling" << G4endl;
    G4double alphamin;
    G4double alphamax;
    //G4double emin;
    G4double emax;
    if (etot > psum[3]) {
      // We should increase the energy
      alphamin = 1;
      //emin = plab[3];
      alphamax = 10;
      emax = 0;
      // Calculate emax
      for (G4int i = 0; i < n; ++i) {
        G4LorentzVector ppp = secs[i].Get4Momentum().boost(-cms);
        G4ThreeVector lll = alphamax * ppp.vect();
        ppp.set(lll,
                std::sqrt(lll.mag2() + secs[i].GetDefinition()->GetPDGMass() * secs[i].GetDefinition()->GetPDGMass()));
        ppp = ppp.boost(cms);
        emax += ppp[3];
      }
    } else {
      // We should decrease the energy
      alphamin = 0;
      //emin = plab[3];
      alphamax = 1;
      emax = psum[3];
    }
    while (std::abs(alphamin - alphamax) > 1e-6) {
      G4double alphamid = 0.5 * (alphamin + alphamax);
      G4double emid = 0;
      for (G4int i = 0; i < n; ++i) {
        G4LorentzVector ppp = secs[i].Get4Momentum().boost(-cms);
        G4ThreeVector lll = alphamid * ppp.vect();
        ppp.set(lll,
                std::sqrt(lll.mag2() + secs[i].GetDefinition()->GetPDGMass() * secs[i].GetDefinition()->GetPDGMass()));
        ppp = ppp.boost(cmsori);
        emid += ppp[3];
      }
      if (etot > emid)
        alphamin = alphamid;
      else
        alphamax = alphamid;
    }
    for (G4int i = 0; i < n; ++i) {
      G4LorentzVector ppp = secs[i].Get4Momentum().boost(-cms);
      G4ThreeVector lll = 0.5 * (alphamin + alphamax) * ppp.vect();
      ppp.set(lll,
              std::sqrt(lll.mag2() + secs[i].GetDefinition()->GetPDGMass() * secs[i].GetDefinition()->GetPDGMass()));
      ppp = ppp.boost(cms);
      secs[i].Set4Momentum(ppp);
    }
  }

  return true;
}
