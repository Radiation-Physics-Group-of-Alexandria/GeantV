// @(#)root/base:$Id: $
// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPXsec
#define ROOT_TPXsec


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPXSec                                                               //
//                                                                      //
// X-section for G5 per particle                                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TDatabasePDG.h"
#include "TPartIndex.h"

class TPXsec {
public:
   TPXsec();
   TPXsec(Int_t pdg, Int_t nxsec);
   virtual ~TPXsec();
   void Print(Option_t *opt="") const;
   const char* Name() const {return TDatabasePDG::Instance()->GetParticle(fPDG)->GetName();}
   Bool_t SetPart(Int_t pdg, Int_t nxsec);
   Bool_t SetPartXS(const Float_t xsec[], const Int_t dict[]);
   Bool_t SetPartIon(const Float_t dedx[]);
   Bool_t SetPartMS(const Float_t angle[], const Float_t ansig[],
		    const Float_t length[], const Float_t lensig[]);
   Int_t PDG() const {return fPDG;}
   Float_t XS(Int_t rindex, Double_t en) const;
   Bool_t XS_v(Int_t npart, Int_t rindex, const Double_t en[], Double_t lam[]) const;
   Float_t DEdx(Double_t en) const;
   Bool_t MS(Double_t en, Float_t &ang, Float_t &asig, 
	     Float_t &len, Float_t &lsig) const;
   Bool_t Resample();
   Bool_t Prune();
   Int_t SampleReac(Double_t en) const;
   Int_t SampleReac(Double_t en, Double_t randn)  const;

   void Dump() const;
   void Interp(Double_t egrid[], Float_t value[], Int_t nbins, 
	       Double_t eildelta, Int_t stride, Double_t en, Float_t result[]);
   
   static void SetVerbose(Int_t verbose) {fVerbose=verbose;}
   static Int_t GetVerbose() {return fVerbose;}
private:
   TPXsec(const TPXsec&); // Not implemented
   TPXsec& operator=(const TPXsec&); // Not implemented

   static Int_t    fVerbose;       // Controls verbosity level

   Int_t           fPDG;           // particle pdg code
   Int_t           fNEbins;        // number of energy bins
   Int_t           fNCbins;        // number of energy bins for dEdx and MS
   Int_t           fNXsec;         // number of reactions
   Int_t           fNTotXs;        // tot size of fTotXs
   Int_t           fNXSecs;        // tot size of fXSecs
   Double_t        fEmin;          // Min energy of the energy grid
   Double_t        fEmax;          // Max energy of the energy grid
   Double_t        fEilDelta;      // logarithmic energy delta
   const Double_t *fEGrid;         //![fNEbins] energy grid
   Float_t        *fMSangle;       // [fNCbins] table of MS average angle
   Float_t        *fMSansig;       // [fNCbins] table of MS sigma angle
   Float_t        *fMSlength;      // [fNCbins] table of MS average lenght correction
   Float_t        *fMSlensig;      // [fNCbins] table of MS sigma lenght correction
   Float_t        *fdEdx;          // [fNCbins] table of dE/dx
   Float_t        *fTotXs;         // [fNTotXs] table of total x-sec
   Float_t        *fXSecs;         // [fNXSecs] table of partial x-sec
   Int_t           fRdict[FNPROC]; // reaction dictionary from reaction number to position
                                  // in the X-sec array
   Int_t           fRmap[FNPROC];  // reaction map, from reaction position in the X-sec
                                  // array to the raction number

   ClassDef(TPXsec,1)  //Particle X-secs
};

#endif
