// This class is to doppler broaden the point ENDF data. It has been derived from F. Carminati's Fortran Function
// Author: Dr. Harphool Kumawat
// Email: harphool@barc.gov.in; harphool.kumawat@cern.ch
// date of creation: March 22, 2016

#include <math.h>
#include <iostream>
#include "TNudyEndfDoppler.h"
#include <iomanip>
#include <algorithm>

using namespace Nudy;
using namespace NudyPhysics;

#ifdef USE_ROOT
ClassImp(TNudyEndfDoppler)
#endif

    TNudyEndfDoppler::TNudyEndfDoppler()
{
}

TNudyEndfDoppler::TNudyEndfDoppler(double isigDiff, double aw, double t1, double t2, std::vector<double> &x1,
                                   std::vector<double> &x2)
    : F0K2P(0), F1K2P(0), F2K2P(0), F3K2P(0), F4K2P(0)
{
#define XNEPSM(S) fabs(ONE - TW3 * sqrt(ONE / 3.0) * pow((ONE + S * (ONE + S)), THH) / (S * (ONE + S)))
#define FTAIL(X, Y)                                                                                                \
  OVSQPI *((1.0 + 2.0 * Y * Y) * sqrt(PI) * HALF * (erf(X - Y) - erf(X + Y)) - (X + Y) * exp(-(X - Y) * (X - Y)) \
            + (X - Y) * exp(-(X + Y) * (X + Y)))
  sigDiff = isigDiff ;
  awri  = aw;
  tk    = t2 - t1;
  if(t1 == t2){
    for (unsigned int j = 0; j < x1.size(); j++) {
      sigma.push_back(x2[j]);  
    } 
  } else {
    ALPHA = awri / ( boltz * tk );
    ncrs  = x1.size();
    while (RATHIG - RATLOW > 1E-7) {
      RATHLF = HALF * (RATLOW + RATHIG);
      HTEST  = XNEPSM(RATHLF);
      if (HTEST < isigDiff) {
	RATLOW = RATHLF;
      } else {
	RATHIG = RATHLF;
      }
    } // end of while loop
    RATMAX = RATLOW * RATLOW;
    EMAX   = x1[ncrs - 1];
    EMIN   = x1[0];
    if ( x1[0] == x1[1] ) {
      x1.erase (x1.begin()+1) ;
      x2.erase (x2.begin()+1) ;
    }
    IPP  = 0;
    size = x1.size();
    jloop = 0 ;
    xss = 0.5 * (2*x2[0]+(x2[1]-x2[0])*(-x1[0])/(x1[1]-x1[0])) ;
    if (xss < 0) xss = 0 ;
    for (int k = 0 ; k < size ; k++) {
//        std::cout <<"Loop begins "<< x1[jloop] <<"  "<< x2[jloop] << std::endl;
      Y2  = x1[jloop] * ALPHA;
      Y   = sqrt(Y2);
      ZKT = sqrt(x1[IPP + 1] * ALPHA) - Y;
      while (ZKT < -ZLIMI) {
	IPP += 1;
	ZKT = sqrt(x1[IPP + 1] * ALPHA) - Y;
	if( ZKT + ZLIMI > 0.05 ){
	  IPP -= 2 ;
	  ZKT = sqrt(x1[IPP + 1] * ALPHA) - Y;
	  break;
	}
      } // end of while loop
      mipp = IPP ;
      XSUM = 0.0;
      FTAIL1 = 0.0;
      KPP  = IPP;
      E2   = x1[KPP];
      S2   = x2[KPP];
      ZK2  = sqrt(E2 * ALPHA) - Y;
      ZK22 = ZK2 * ZK2;
      EXPA = exp(-ZK22);
      F0K2 = erf(ZK2);
      F1K2 = OVSQPI * (ONE - EXPA);
      F2K2 = HALF * F0K2 - OVSQPI * ZK2 * EXPA;
      F3K2 = OVSQPI * (ONE - (1 + ZK22) * EXPA);
      F4K2 = HALF * THH * F0K2 - OVSQPI * ZK2 * (THH + ZK22) * EXPA;
      
      if (Y < ZLIMI) {
	ZK2P  = ZK2 + 2 * Y;
	ZK22P = ZK2P * ZK2P;
	EXPAP = exp(-ZK22P);
	F0K2P = erf(ZK2P);
	F1K2P = OVSQPI * (ONE - EXPAP);
	F2K2P = HALF * F0K2P - OVSQPI * ZK2P * EXPAP;
	F3K2P = OVSQPI * (ONE - (1 + ZK22P) * EXPAP);
	F4K2P = HALF * THH * F0K2P - OVSQPI * ZK2P * (THH + ZK22P) * EXPAP;
      }
      
      while (ZK2 < ZLIMI && KPP < (int) x1.size() - 1) {
	E1  = E2;
	S1  = S2;
	KPP = KPP + 1;
	E2  = x1[KPP];
	S2  = x2[KPP];
	while ( E2 == E1 ) {
	  KPP = KPP + 1;
	  E2  = x1[KPP];
	  S2  = x2[KPP];
	}
	if (E2 - E1 == 0.0 || E2 == 0 || ALPHA == 0.0) {
	  std::cout << "Doppler fails between " << E1 << "  " << E2 << "  " << Y2 << "  " << ALPHA << "  " << awri << std::endl;
	  continue;
	}
	F0K1 = F0K2;
	F1K1 = F1K2;
	F2K1 = F2K2;
	F3K1 = F3K2;
	F4K1 = F4K2;
	ZK2  = sqrt(E2 * ALPHA) - Y;
	ZK22 = ZK2 * ZK2;
	EXPA = exp(-ZK22);
	F0K2 = erf(ZK2);
	F1K2 = OVSQPI * (ONE - EXPA);
	F2K2 = HALF * F0K2 - OVSQPI * ZK2 * EXPA;
	F3K2 = OVSQPI * (ONE - (1 + ZK22) * EXPA);
	F4K2 = HALF * THH * F0K2 - OVSQPI * ZK2 * (THH + ZK22) * EXPA;
	FACT = ONE / (E2 - E1);
	AK   = (E2 * S1 - E1 * S2) * FACT;
	CK   = (S2 - S1) * FACT / ALPHA;
	CKY  = CK * Y;
	CKY2 = CK * Y2;
	XSUM = XSUM + CK * (F4K2 - F4K1) + 4 * CKY * (F3K2 - F3K1) + (AK + 6 * CKY2) * (F2K2 - F2K1) +
	      2 * Y * (AK + 2 * CKY2) * (F1K2 - F1K1) + Y2 * (AK + CKY2) * (F0K2 - F0K1);
	      
	if (Y < ZLIMI) {
	  ZK1P  = ZK2P;
	  F0K1P = F0K2P;
	  F1K1P = F1K2P;
	  F2K1P = F2K2P;
	  F3K1P = F3K2P;
	  F4K1P = F4K2P;
	  ZK2P  = ZK2 + 2 * Y;
	  ZK22P = ZK2P * ZK2P;
	  EXPAP = exp(-ZK22P);
	  F0K2P = erf(ZK2P);
	  F1K2P = OVSQPI * (ONE - EXPAP);
	  F2K2P = HALF * F0K2P - OVSQPI * ZK2P * EXPAP;
	  F3K2P = OVSQPI * (ONE - (1 + ZK22P) * EXPAP);
	  F4K2P = HALF * THH * F0K2P - OVSQPI * ZK2P * (THH + ZK22P) * EXPAP;
	  XSUM  = XSUM - (CK * (F4K2P - F4K1P) - 4 * CKY * (F3K2P - F3K1P) + (AK + 6 * CKY2) * (F2K2P - F2K1P) -
			2 * Y * (AK + 2 * CKY2) * (F1K2P - F1K1P) + Y2 * (AK + CKY2) * (F0K2P - F0K1P));
	} // end of if
      }   // end of if
      // while loop
      if (XSUM < 0) XSUM = 0 ;
//      std::cout << HALF * XSUM / Y2 <<" xsum1 \t"<< xss  <<" FTAIL1 \t"<< 2*x2[0] <<" FTAIL2 \t"<< x2[ jloop ] <<std::endl;
      FTAIL1 = xss * (FTAIL(sqrt(x1[0] * ALPHA), Y) - FTAIL(ZERO, Y));
      XSUM   = XSUM + FTAIL1;
      FTAIL2 = x2[x1.size() - 1] * (FTAIL(sqrt((x1[x1.size() -1] + 0.1 * x1[x1.size() -1])* ALPHA), Y) - FTAIL(sqrt(x1[x1.size() - 1] * ALPHA), Y));
      XSUM   = XSUM + FTAIL2;
      sigma.push_back(HALF * XSUM / Y2);
//        std::cout << HALF * XSUM / Y2 <<" xsum2 \t"<< HALF * FTAIL1 / Y2 <<" FTAIL1 \t"<< HALF * FTAIL2 / Y2 <<" FTAIL2 \t"<< x2[ jloop ] <<std::endl;
      if ( jloop > 0 && k < size - 1) {
	mloop = 0 ;
	recursionLinear1(x1, x2, x1[ jloop - 1 ], x2[ jloop - 1 ], sigma[ jloop - 1 ], x1[ jloop ], x2[ jloop ], sigma[ jloop ]);
	jloop += mloop ;
      }
      jloop++ ;	
    } // end of for loop
  }
}

double TNudyEndfDoppler::recursionLinear1(std::vector<double> &x1,
					  std::vector<double> &x2,
					  double x,
					  double y,
					  double sig,
					  double xd,
					  double yd,
					  double sigd ) {
#define XNEPSM(S) fabs(ONE - TW3 * sqrt(ONE / 3.0) * pow((ONE + S * (ONE + S)), THH) / (S * (ONE + S)))
#define FTAILX(X, Y)                                                                                                \
  OVSQPI *((1.0 + 2.0 * Y * Y) * sqrt(PI) * HALF * (erf(X - Y) - erf(X + Y)) - (X + Y) * exp(-(X - Y) * (X - Y)) + \
           (X - Y) * exp(-(X + Y) * (X + Y)))
   if ( y <= 0.0 && yd <= 0.0) return 0;
  double mid 		 = 0.5 * ( x + xd );
  double sigmid1         = y + (yd - y) * (mid - x) / (xd - x);
  double sigmid2         = sig + (sigd - sig) * (mid - x) / (xd - x);
  std::vector<double>::iterator itx;
  itx = std::find (x1.begin(), x1.end(), x);
  int xindex = itx - x1.begin() ;
  x1.insert(x1.begin() + xindex + 1, 1 , mid);
  x2.insert(x2.begin() + xindex + 1 , 1 , sigmid1);
  double sigmid3 = broadMore ( x1, x2, mid ) ;
  double errmid = fabs(sigmid2 - sigmid3)/sigmid3 ;
//  if ( errmid >= sigDiff ) {
  if (fabs(sig/sigd -1) > 1E-2 || errmid >= sigDiff ) {
//  if (m1 > 0 && m3 > 0 && mdiff > 1E-3  && mdiff != 1) {
    sigma.insert (sigma.begin() + xindex + 1, 1 , sigmid3 );
    mloop++ ;
//      std::cout << x <<"  \t"<< xd <<"  \t"<< errmid <<"  \t"<<sigmid3<<"  \t"<<sigmid2 << std::endl;
//      std::cout << sig <<"  \t"<< sigd <<"  \t"<< mid << std::endl;
//      std::cout << x <<"  \t"<< y <<"  \t"<< xd <<"  \t"<< yd << std::endl;
    recursionLinear1( x1, x2, x , y, sig, mid, sigmid1, sigmid3 ) ; 
    recursionLinear1( x1, x2, mid, sigmid1, sigmid3, xd, yd, sigd ) ; 
    return 0;
  } else {
    x1.erase(x1.begin() + xindex + 1);
    x2.erase(x2.begin() + xindex + 1);
    return 0;
  }  
  return 0;
}
double TNudyEndfDoppler::broadMore (std::vector<double> &x1 , 
				    std::vector<double> &x2 ,
				    double ixp ) {
  
  IPP = mipp ;
  KPP = IPP ;
  ncrs  = x1.size() ;
  Y2  = ixp * ALPHA ;
  Y   = sqrt(Y2) ;
  ZKT = sqrt(x1[IPP + 1] * ALPHA) - Y ;
  while (ZKT > -ZLIMI) {
    IPP -= 1;
    ZKT = sqrt( x1 [ IPP - 1 ] * ALPHA ) - Y;
    if(IPP <= 0)
    {
      IPP = 0;
      break;
    }
  } // end of while loop
  XSUM = 0.0;
  FTAIL1 = 0.0;
  KPP  = IPP;
  E2   = x1 [ KPP ];
  S2   = x2 [ KPP ];
  ZK2  = sqrt(E2 * ALPHA) - Y;
  ZK22 = ZK2 * ZK2;
  EXPA = exp(-ZK22);
  F0K2 = erf(ZK2);
  F1K2 = OVSQPI * (ONE - EXPA);
  F2K2 = HALF * F0K2 - OVSQPI * ZK2 * EXPA;
  F3K2 = OVSQPI * (ONE - (1 + ZK22) * EXPA);
  F4K2 = HALF * THH * F0K2 - OVSQPI * ZK2 * (THH + ZK22) * EXPA;
  if (Y < ZLIMI) {
    ZK2P  = ZK2 + 2 * Y;
    ZK22P = ZK2P * ZK2P;
    EXPAP = exp(-ZK22P);
    F0K2P = erf(ZK2P);
    F1K2P = OVSQPI * (ONE - EXPAP);
    F2K2P = HALF * F0K2P - OVSQPI * ZK2P * EXPAP;
    F3K2P = OVSQPI * (ONE - (1 + ZK22P) * EXPAP);
    F4K2P = HALF * THH * F0K2P - OVSQPI * ZK2P * (THH + ZK22P) * EXPAP;
  }
  while (ZK2 < ZLIMI && KPP <  (int)x1.size() - 1) {
    E1  = E2;
    S1  = S2;
    KPP = KPP + 1;
    E2  = x1[KPP];
    S2  = x2[KPP];
    while ( E2 == E1 ) {
      KPP = KPP + 1;
      E2  = x1[KPP];
      S2  = x2[KPP];
    }
    if (E2 - E1 == 0.0 || Y2 == 0.0 || ALPHA == 0.0) {
      std::cout << "Doppler fails " << E1 << "  " << E2 << "  " << Y2 << "  " << ALPHA << "  " << awri << std::endl;
      continue;
    }
    F0K1 = F0K2;
    F1K1 = F1K2;
    F2K1 = F2K2;
    F3K1 = F3K2;
    F4K1 = F4K2;
    ZK2  = sqrt(E2 * ALPHA) - Y;
    ZK22 = ZK2 * ZK2;
    EXPA = exp(-ZK22);
    F0K2 = erf(ZK2);
    F1K2 = OVSQPI * (ONE - EXPA);
    F2K2 = HALF * F0K2 - OVSQPI * ZK2 * EXPA;
    F3K2 = OVSQPI * (ONE - (1 + ZK22) * EXPA);
    F4K2 = HALF * THH * F0K2 - OVSQPI * ZK2 * (THH + ZK22) * EXPA;
    FACT = ONE / (E2 - E1);
    AK   = (E2 * S1 - E1 * S2) * FACT;
    CK   = (S2 - S1) * FACT / ALPHA;
    CKY  = CK * Y;
    CKY2 = CK * Y2;
    XSUM = XSUM + CK * (F4K2 - F4K1) + 4 * CKY * (F3K2 - F3K1) + (AK + 6 * CKY2) * (F2K2 - F2K1) +
	  2 * Y * (AK + 2 * CKY2) * (F1K2 - F1K1) + Y2 * (AK + CKY2) * (F0K2 - F0K1);
    if (Y < ZLIMI) {
      ZK1P  = ZK2P;
      F0K1P = F0K2P;
      F1K1P = F1K2P;
      F2K1P = F2K2P;
      F3K1P = F3K2P;
      F4K1P = F4K2P;
      ZK2P  = ZK2 + 2 * Y;
      ZK22P = ZK2P * ZK2P;
      EXPAP = exp(-ZK22P);
      F0K2P = erf(ZK2P);
      F1K2P = OVSQPI * (ONE - EXPAP);
      F2K2P = HALF * F0K2P - OVSQPI * ZK2P * EXPAP;
      F3K2P = OVSQPI * (ONE - (1 + ZK22P) * EXPAP);
      F4K2P = HALF * THH * F0K2P - OVSQPI * ZK2P * (THH + ZK22P) * EXPAP;
      XSUM  = XSUM - (CK * (F4K2P - F4K1P) - 4 * CKY * (F3K2P - F3K1P) + (AK + 6 * CKY2) * (F2K2P - F2K1P) -
		    2 * Y * (AK + 2 * CKY2) * (F1K2P - F1K1P) + Y2 * (AK + CKY2) * (F0K2P - F0K1P));
    } // end of if
  }   // end of if
  // while loop
  if (XSUM < 0) XSUM = 0 ;
  FTAIL1 = xss * (FTAILX(sqrt(x1[0] * ALPHA), Y) - FTAILX(ZERO, Y));
  XSUM   = XSUM + FTAIL1;
  FTAIL2 = x2[x1.size() - 1] * (FTAILX(sqrt((x1[x1.size() -1] + 0.1 * x1[x1.size() -1])* ALPHA), Y) - FTAILX(sqrt(x1[x1.size() - 1] * ALPHA), Y));
  XSUM   = XSUM + FTAIL2;
  return HALF * XSUM / Y2;
}
