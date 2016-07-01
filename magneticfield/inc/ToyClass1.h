#ifndef _TOYCLASS1_H_
#define _TOYCLASS1_H_

#include <iostream>
#include "AlignedBase.h"
#include <Vc/Vc>
#include <stdlib.h>
#include "base/Vector3D.h" //To utlimately include kVc, kScalar backends etc


using namespace std;
using Backend = vecgeom::kVc;

class ToyClass1 : public AlignedBase
{
  typedef typename vecgeom::kVc::bool_v      Bool_v;
  typedef typename vecgeom::kVc::precision_v Double_v;
  typedef typename vecgeom::kVc::int_v       Int_v;
public:
  ToyClass1(){};

  //used to set total number of tracks user might give
  //Can be modified later. Makes things easier for now.
  ToyClass1(int n);

  ~ToyClass1();

  //Main method here
  void ToyMethod(double nTracks[], double finalResults[]);

  //Called by ToyMethod; kinda like OneGoodStep
  void SimpleMethod(Double_v &fPreProcLane, 
                    Double_v &outputSimpleMethod);

  //Called in constructor. Set input no. of tracks
  void SetInputTracks(int n);

  //Initialize the Vc vector in beginning if not initialized. Assume 
  //uninitialized here. Later on replacement when work done for it 
  void InitializeVcVector(double nTracks[]);

  void PreProcess(double nTracks[], 
                  double fPreprocessedData[]);

  Bool_v IsGoodStep(Double_v outputSimpleMethod);

  //Print statements for debugging
  void PrintCurrentData(Double_v fPreProcLane, 
                        Double_v outputSimpleMethod,
                        Bool_v Done                );
  void PrintCurrentData(Double_v outputSimpleMethod,
                        Bool_v Done                );

private:
  int kSizeOfVcVector = 4; //can use vecgeom::kVectorSize
  int fInputTotalTracks;
  double kStep = 0.1;

  //Stores indices, used to take care of sth sth.. random drops and pickups
  // Int_v fIndex;
  int *fIndex;

  //Stores scalar preprocessed data 
  //probably needs to be vectorized but leave that for later as of now 
  //this is kinda mix of scalar and vector version
  //but let it be for the time being 
  double *fPreprocessedData;

  Double_v fNoGoodStepsLane;
  Double_v fYVectorLane;
  Double_v fPreProcLane;
  Double_v fStepStateLane;

};


ToyClass1::ToyClass1(int n)
{
  SetInputTracks(n);
  fPreprocessedData = new double[n];
  fIndex       = new int[kSizeOfVcVector];
}

ToyClass1::~ToyClass1()
{
  delete fPreprocessedData;
  cout<<"----ToyClass1 destructor being called"<<endl;
}

void ToyClass1::SetInputTracks(int n)
{
  fInputTotalTracks = n;
}

void ToyClass1::InitializeVcVector(double nTracks[])
{
  for (int i = 0; i < kSizeOfVcVector; ++i)
  {
    fPreProcLane    [i] = nTracks[i];
    fIndex          [i] = i;    
    fStepStateLane  [i] = kStep; 
    fNoGoodStepsLane[i] = 0;
  }
}

void ToyClass1::PreProcess(double nTracks[], double fPreprocessedData[])
{  
  double factor = 0.4;

  for (int i = 0; i < fInputTotalTracks; ++i)
  {
    fPreprocessedData[i] = nTracks[i] * factor;
  }
}

void ToyClass1::PrintCurrentData(typename vecgeom::kVc::precision_v fPreProcLane, 
                                 typename vecgeom::kVc::precision_v outputSimpleMethod, 
                                 typename vecgeom::kVc::bool_v      Done                )
{
  cout<<"\n----Currently:"<<endl;
  cout<<"----fPreProcLane is      : "<<fPreProcLane      <<endl;
  cout<<"----outputSimpleMethod is: "<<outputSimpleMethod<<endl;
  cout<<"----Done is              : "<<Done              <<endl;
}

void ToyClass1::PrintCurrentData(typename vecgeom::kVc::precision_v outputSimpleMethod, 
                                 typename vecgeom::kVc::bool_v      Done               )
{
  cout<<"\n----Currently:"<<endl;
  cout<<"----outputSimpleMethod is: "<<outputSimpleMethod<<endl;
  cout<<"----Done is              : "<<Done<<endl;
}


void ToyClass1::SimpleMethod(typename vecgeom::kVc::precision_v &fPreProcLane, 
                             typename vecgeom::kVc::precision_v &outputSimpleMethod)
{
  // srand(time(NULL));
  for (int i = 0; i < kSizeOfVcVector; ++i)
  {
    outputSimpleMethod[i] =  (float) rand()/(RAND_MAX) ;
  }
  // outputSimpleMethod = fPreProcLane/2.; //some random processing
  fPreProcLane /= 3.;

}

typename Backend::bool_v 
ToyClass1::IsGoodStep(typename vecgeom::kVc::precision_v outputSimpleMethod)
{
  Bool_v isGoodStep = (outputSimpleMethod > fStepStateLane - kStep ) && (outputSimpleMethod < fStepStateLane);
  return isGoodStep;
}

void ToyClass1::ToyMethod(double nTracks[], double finalResults[])
{
  PreProcess(nTracks, fPreprocessedData);

  InitializeVcVector(fPreprocessedData);

  int trackNextInput = 4; 

  //now start working from preprocessed data
  typedef typename vecgeom::kVc::bool_v Bool_v;
  typedef typename vecgeom::kVc::precision_v Double_v;
  Bool_v isDone(0.), isGoodStep(0.);
  Double_v outputSimpleMethod;


  while(!(vecgeom::IsFull(isDone)) || (trackNextInput < fInputTotalTracks))
  {

    SimpleMethod(fPreProcLane, outputSimpleMethod);
    //Check if between fstate-0.1 and fstate
    //fstate beginning from 0.1
    //So first step is between 0 and 0.1
    // isGoodStep = (outputSimpleMethod > fStepStateLane - kStep ) && (outputSimpleMethod < fStepStateLane);
    //Function IsGoodStep instead of checking here
    isGoodStep = IsGoodStep(outputSimpleMethod);

    //if it is a good step, increment state by 0.1 and increase no. of good steps
    vecgeom::MaskedAssign(isGoodStep && !isDone, fStepStateLane + kStep, &fStepStateLane  );
    vecgeom::MaskedAssign(isGoodStep, fNoGoodStepsLane + 1,   &fNoGoodStepsLane);

    //now, if no. of good steps =9/10, then isDone = true and insert a new track and say isDone = false 
    //for the new track
    // if no track left, then leave isDone = true;

    isDone = (fNoGoodStepsLane >=10);

    PrintCurrentData(outputSimpleMethod, isDone);
    cout<<"----isGoodStep is      : "<<isGoodStep      <<endl;
    cout<<"----fStepStateLane is  : "<<fStepStateLane  <<endl;
    cout<<"----fNoGoodStepsLane is: "<<fNoGoodStepsLane<<endl;

    //what if all are not done?
    //if none is done, then simply continue
    //But if say one is done, then we need to insert new tracks in it's position
    //This executes only if at least one is done 
    //Do it only if anything is left in nTracks, otherwise no extra tracks left
    // cout<<  !(vecgeom::IsFull(isDone)) <<" "<< vecgeom::IsEmpty(isDone) <<endl;
    // if( !(vecgeom::IsFull(isDone)) && !(vecgeom::IsEmpty(isDone)))
    if(!(vecgeom::IsEmpty(isDone)))
    {
      cout<<"here1"<<endl;
      for (int i = 0; i < kSizeOfVcVector; ++i)
      {
        cout<<"here2"<<endl;
        if(isDone[i]==1 && fIndex[i] != -1) 
        {
          cout<<"here3"<<endl;
          finalResults[fIndex[i]] = outputSimpleMethod[i]; //store the output

          if(trackNextInput < fInputTotalTracks)
          { 
            cout<<"\n----trackNextInput is: "<<trackNextInput<<endl;
            //insert new track
            //And store output of this one in finalResults (in scalar array for now)
            //Probably need to define an index to take care of these random pickups and drops
            fPreProcLane    [i]  = nTracks[trackNextInput]; //sending in next one
            fIndex          [i]  = trackNextInput;
            fStepStateLane  [i]  = kStep;
            fNoGoodStepsLane[i]  = 0;
            isDone          [i]  = 0;
            trackNextInput++;
          }
          else
          {
            cout<<"here4"<<endl;
            isDone [i] =  1;
            fIndex [i] = -1;
            //need to stop processing for other things or do something else of the kind
          }
        }
      }
    }

  }

  cout<<"----Exited while loop"<<endl;

  cout<<"----here we are at the end"<<endl;
  cout<<"\n----Output obtained is: "<<endl;
  for (int i = 0; i < fInputTotalTracks; ++i)
  {
    cout<<i<<" "<<finalResults[i]<<" "<<endl;
  }

  return ;
}





#endif