// @(#)root/base:$Id: $
// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPartIndex
#define ROOT_TPartIndex

#define DICLEN 12

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPartIndex                                                           //
//                                                                      //
// Particle index singleton for various particle translation functions  //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <RTypes.h>

class TPartIndex {

public:
   const TPartIndex* I() {if(!fgPartIndex) fgPartIndex=this; return fgPartIndex;}
   virtual ~TPartIndex();

private:
   TPartIndex(): fNpart(0), fNpartReac(0), fPDG(0), fPDGReac(0) {}
   static TPartIndex  *fgPartIndex;
   Int_t          fNpart; // Total number of particles
   Int_t          fNpartReac; // Total number of particles with reactions
   Int_t          fNReac; //
   Short_t       *fPDG; // [fNpart] Correspondence PDG <-> particle number (from G4)
   Short_t       *fPDGReac; // [fNpartReac] Correspondence PDG <-> particle number with reac
   Char_t        *fRdict[12]; // [fNReac]

   ClassDef(TPartIndex,1)  // Particle Index

};


#endif
