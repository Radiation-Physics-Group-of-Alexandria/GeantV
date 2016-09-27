#include "GUAliasSampler.h"

#include "MaterialHandler.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECCORE_ATT_HOST
GUAliasSampler::GUAliasSampler(Random_t *states, int threadId, double incomingMin, double incomingMax,
                               int numEntriesIncoming, // for 'energy' (or log) of projectile
                               int numEntriesSampled)
    : fRandomState(states), fThreadId(threadId), fIncomingMin(incomingMin), fIncomingMax(incomingMax),
      fInNumEntries(numEntriesIncoming), fLogIncomingMin(math::Log(incomingMin)),
      fInverseBinIncoming(numEntriesIncoming / (incomingMax - incomingMin)),
      fInverseLogBinIncoming(numEntriesIncoming / (math::Log(incomingMax) - fLogIncomingMin)),
      fSampledNumEntries(numEntriesSampled), fInverseBinSampled(1.0 / numEntriesSampled)
{
  int nelements = MaterialHandler::Instance()->GetNumberOfElements();
  int ngrid = (fInNumEntries + 1) * fSampledNumEntries;
  fAliasTableManager = new GUAliasTableManager(nelements, ngrid);
}

VECCORE_ATT_HOST_DEVICE
GUAliasSampler::GUAliasSampler(Random_t *states, int threadId, double incomingMin, double incomingMax,
                               int numEntriesIncoming, // for 'energy' (or log) of projectile
                               int numEntriesSampled, GUAliasTableManager *tableManager)
    : fRandomState(states), fThreadId(threadId), fIncomingMin(incomingMin), fIncomingMax(incomingMax),
      fInNumEntries(numEntriesIncoming), fLogIncomingMin(math::Log(incomingMin)),
      fInverseBinIncoming(numEntriesIncoming / (incomingMax - incomingMin)),
      fInverseLogBinIncoming(numEntriesIncoming / (math::Log(incomingMax) - fLogIncomingMin)),
      fSampledNumEntries(numEntriesSampled),
      fInverseBinSampled(1.0 / numEntriesSampled) // Careful - convention build / use table!
{
  fAliasTableManager = tableManager;
}

VECCORE_ATT_HOST_DEVICE
GUAliasSampler::~GUAliasSampler()
{
#if !defined(VECCORE_CUDA)
  if (fAliasTableManager)
    delete fAliasTableManager;
#endif
}

VECCORE_ATT_HOST_DEVICE
void GUAliasSampler::PrintTable()
{
  printf("Incoming Min= %g , Max= %g , numEntries= %d \n", fIncomingMin, fIncomingMax, fInNumEntries);

  if (fAliasTableManager->GetNumberOfElements() > 0) {
    for (int i = 0; i < fAliasTableManager->GetNumberOfElements(); ++i) {
      GUAliasTable *tb = fAliasTableManager->GetAliasTable(i);
      printf("GUAliasSampler fAliasTable = %p fNGrid = %d value=(%d,%f,%f)\n", tb, tb->fNGrid, tb->fAlias[1],
             tb->fProbQ[1], tb->fpdf[1]);
    }
  }
  else {
    printf("GUAliasSampler fAliasTableManager is empty\n");
  }
}

VECCORE_ATT_HOST
void GUAliasSampler::BuildAliasTable(int Zelement, const double *pdf)
{
  // Build alias and alias probability
  //
  // Reference: (1) A.J. Walker, "An Efficient Method for Generating Discrete
  // Random Variables with General Distributions" ACM Trans. Math. Software, 3,
  // 3, 253-256 (1977) (2) A.L. Edwards, J.A. Rathkopf, and R.K. Smidt,
  // "Extending the Alias Monte Carlo Sampling Method to General Distributions"
  // UCRL-JC-104791 (1991)
  //
  // input : fInNumEntries       multiplicity of alias tables)
  //         fSampledNumEntries (dimension of discrete outcomes)
  //         pdf[fInNumEntries x fSampledNumEntries] (probability density function)
  // output: a[fInNumEntries x fSampledNumEntries]   (alias)                        ----> :  a[fSampledNumEntries]
  //         q[fInNumEntries x fSampledNumEntries]   (non-alias probability)        ----> :  q[fSampledNumEntries]
  //

  // temporary array
    
  int *a = (int *)malloc(fSampledNumEntries * sizeof(int));             // : This variable is not used
  double *ap = (double *)malloc(fSampledNumEntries * sizeof(double));   // : what's this? why not q?

  // likelihood per equal probable event
  const double cp = 1.0 / fSampledNumEntries;
  GUAliasTable *table = new GUAliasTable((fInNumEntries + 1) * fSampledNumEntries);

#ifdef ALIAS_DEBUG
    std::cout<<"\n***** GUAliasSampler::BuildAliasTable debug START.\n";
    std::cout<<"Alias Table created ["<<fInNumEntries + 1<<"]x["<<fSampledNumEntries<<"]\n";
    std::cout<<"Likelihood per equal probable events: "<<cp<<"\n";
    int *selectedAsDonor=(int *)malloc(fSampledNumEntries * sizeof(int));
    int *selectedAsRecipient=(int *)malloc(fSampledNumEntries * sizeof(int));
    int *myDonors= (int *)malloc(fSampledNumEntries*fSampledNumEntries * sizeof(int));
    int *myRecipients= (int *)malloc(fSampledNumEntries*fSampledNumEntries * sizeof(int));
    for(int i=0; i< fSampledNumEntries; i++)
    {
        selectedAsDonor[i]=0;
        selectedAsRecipient[i]=0;
        for(int j=0; j< fSampledNumEntries; j++)
        {
            myDonors[i*fSampledNumEntries+j]=-1;
            myRecipients[i*fSampledNumEntries+j]=-1;
        }
    }
#endif
    /*
    ERRORS and CONSIDERATIONS:
     - when I check for a donor, it's not said that simultaneusly I find a recipient and viceversa. This leads to errors. In fact, due to rounding errors and errors due to the use of doubles, one of the two lists could be empty, while the other is not.
     - when a bin has probability=1 then ki must be initialized to i, ki=i, sensitive info.
     - when I control the probabilities I have to keep into consideration a small epsilon (error) due to floating point operations
     - operations of <= don't make sense, see point below
     - need to understand why the dimension of the grid is [n+1]*[m] and not [n]*[m]
     - it's probably necessary to scale the pdf(s) at the beginning in order to be compared with the equal likelihood. Standard implementation of the alias table do this at the beginning to avoid errors.
     - the probability of the last/lasts bin/s, that will be selected as a donor, must be set to 1.
     */
  
for (int ir = 0; ir <= fInNumEntries; ++ir) { //check this <=
    // copy and initialize
    for (int i = 0; i < fSampledNumEntries; ++i) {

      int ipos = ir * fSampledNumEntries + i;
      table->fpdf[ipos] = pdf[ipos]; // Copy of the original pdf

      a[i] = i; // this variable is not used
      ap[i] = pdf[ipos];
    }

    // O(n) iterations
    int iter = fSampledNumEntries; // Iterate for all the bins

#ifdef ALIAS_DEBUG
      double sum=0;
      for(int index=0; index<fSampledNumEntries; ++index)
      {
          sum=sum+ap[index];
          std::cout<<"Starting n-AliasProb["<<index<<"]="<<ap[index]<<"\n";
      }
      std::cout<<"Equal likelihood value: "<<cp<<" and sum of all the pdf: "<<sum<<"\n";
#endif

    do {
#ifdef ALIAS_DEBUG
        std::cout<<"\nIteration START: "<<iter<<"\n";
        for (int j = 0; j < fSampledNumEntries; ++j) {
                std::cout<<"n-AliasProb["<<j<<"]="<<ap[j]<<"\n";
            }
#endif
      int donor = 0;
      int recip = 0;

      // A very simple search algorithm
      for (int j = donor; j < fSampledNumEntries; ++j) {
        if (ap[j] >= cp) { //If the pdf is > or = I became a donor! We are comparing double, the equal doesn't make sense. And also if the pdf is equal to cp, the item shouldn't be selected as a donor.
          //if (nap[j] > cp) { //correction
          donor = j;
#ifdef ALIAS_DEBUG
            selectedAsDonor[j]++;
            std::cout<<"bin : "<<j<<" selected as a donor\n";
#endif
            break;
        }
      }

      for (int j = recip; j < fSampledNumEntries; ++j) { //why  recip and not =0?
        if (ap[j] > 0.0 && ap[j] < cp) {
          recip = j;
#ifdef ALIAS_DEBUG
            selectedAsRecipient[j]++;
            if(selectedAsRecipient[j] > 1){
                //std::cout<<"Error! Bin "<<j<<", selected for the second time as recipient. Exiting.\n";
                //return;
            }
            std::cout<<"bin : "<<j<<", selected as a recipient\n";
#endif

          break;
        }
      }
    //PROBLEM: if a bin is selected as a donor and there are no recipient left, the recipient becomes automatically the bin zero. this is not exactly what we want, because bin zero could have been selected before as a recipient, then we would lose some info.
        // alias index and non-alias probability

      table->fAlias[ir * fSampledNumEntries + recip] = donor; // l'alias di recip é il donor
      table->fProbQ[ir * fSampledNumEntries + recip] = fSampledNumEntries * ap[recip]; 
    
 
#ifdef ALIAS_DEBUG
        std::cout<<"Alias index: table->fAlias["<<ir * fSampledNumEntries + recip<<"]: "<<table->fAlias[ir * fSampledNumEntries + recip]<<"\n";
        std::cout<<"Non-alias probability: table->fProbQ["<<ir * fSampledNumEntries + recip<<"]: "<<table->fProbQ[ir * fSampledNumEntries + recip]<<"\n";
        std::cout<<"Donor:  ap["<<donor<<"] before= "<<ap[donor]<<"\n";
        std::cout<<"Recip:  ap["<<recip<<"] before= "<<ap[recip]<<"\n";
#endif

      // update pdf
      ap[donor] = ap[donor] - (cp - ap[recip]);// subtract the portion that was missing in recip to reach cp
      ap[recip] = 0.0;// in this way this bin will never be selected again as a recip
      --iter;
#ifdef ALIAS_DEBUG
        std::cout<<"End Iteration: "<<iter+1<<"\n";
        std::cout<<"Donor: "<<donor<<" ap[donor] after= "<<ap[donor]<<"\n";
        std::cout<<"Recip: "<<recip<<" ap[recip] after= "<<ap[recip]<<"\n";
#endif

    } while (iter > 0);
  }
    
    

#ifdef ALIAS_DEBUG
    std::cout<<"\n** END ITERATIONS ** \n\n";

    for(int index=0; index<fSampledNumEntries; ++index)
    {
        std::cout<<"n-AliasProb["<<index<<"]="<<ap[index]<<"\n";
        std::cout<<"selectedAsDonor["<<index<<"]="<<selectedAsDonor[index]<<" times\n";
        std::cout<<"selectedAsRecipient["<<index<<"]="<<selectedAsRecipient[index]<<" times\n";
        
    }
    //At the end we expect that all the probabilities in the table are positive or equal to zero. This check can be added.
    
    for (int ir = 0; ir <= fInNumEntries; ++ir) {
        for (int i = 0; i < fSampledNumEntries; ++i) {
            std::cout<<" table->fAlias["<<ir * fSampledNumEntries + i<<"]= "<< table->fAlias[ir * fSampledNumEntries + i]<<"\n";
            std::cout<<" table->fProbQ["<<ir * fSampledNumEntries + i<<"]= "<< table->fProbQ[ir * fSampledNumEntries + i]<<"\n";
        }
    }
    
#endif
  fAliasTableManager->AddAliasTable(Zelement, table);

  free(a);
  free(ap);
}

} // end namespace impl
} // end namespace vecphys
