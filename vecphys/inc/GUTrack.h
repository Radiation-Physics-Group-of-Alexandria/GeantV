#include "base/SystemOfUnits.h"

#ifndef GUTRACK_H
#define GUTRACK_H 1

struct GUTrack {
  int status;
  int particleType;
  int id;       // counter
  int parentId; // id of parent
  int proc;     // index of physics process
  double x;     // x position - rarely relevant
  double y;     // y position - ditto
  double z;     // z position - ditto
  double px;
  double py;
  double pz;
  double E;
  double q;      // charge ?
  double nint;   // number of interaction length left
  double lambda; // interaction length
  double s;      // step length ??

 public:
  GUTrack() : status(-1), particleType(0), id(-1), parentId(-2), proc(-1), x(0.), y(0.), z(0.),
     px(0.), py(0.), pz(0.), E(0.0), q(0.0), nint(-1), s(0.0) {}
   
  // static const int  printMode= 1;
  // Streaming operator   
  friend std::ostream& operator<< ( std::ostream& os, const GUTrack& t )
  { return t.StreamInfo(os); }   

  std::ostream& StreamInfo(std::ostream& os) const;
};

// Stream object contents to an output stream

inline
std::ostream& GUTrack::StreamInfo(std::ostream& os) const
{
  const int flpPrec= 4;
  // using CLHEP::cm;  
  using CLHEP::keV;
  using CLHEP::MeV;
  using CLHEP::GeV;  
  
  int oldPrec = os.precision(flpPrec); // For charge
  // os << " guTrack: " // "-----------------------------------------------------------\n"
     // << " *** Dump of GUTrack - " << " ***\n"
     // << " ==========================================================\n";
  os << " id= " << id
     << " type: " << particleType
     << " q= " <<  int(q)
     << " parent = " << parentId
     << " proc chosen= " << proc << "\n";
  // int oldPrec = os.precision(flpPrec);  
  /* os << " Position: (in cm)      : "
     << "    X: " << x/cm << " "
     << "    Y: " << y/cm << " "
     << "    Z: " << z/cm << " \n"
   */
  if( E > 1.0 * GeV ) {
     os << " Momentum: ( in GeV/c ) : "
        << "   pX: " << px / GeV << " "
        << "   pY: " << py / GeV << " "
        << "   pZ: " << pz / GeV << " "
        << " Energy = " << E / GeV << " GeV\n";     
  }
  else if( E > 1.0 * MeV ) {
     os << " Momentum: ( in MeV/c ) : "
        << "   pX: " << px / MeV << " "
        << "   pY: " << py / MeV << " "
        << "   pZ: " << pz / MeV << " "
        << " Energy = " << E / MeV << " MeV\n";
  }
  else{
     os << " Momentum: ( in keV/c ) : "
        << "   pX: " << px / keV << " "
        << "   pY: " << py / keV << " "
        << "   pZ: " << pz / keV << " "
        << " Energy = " << E / keV << " KeV\n";     
  }
  /* os << "  # int-len left = " << nint
     << "  lambda = " << lambda
     << "  s = " << s << "\n" // std::endl
     << "-----------------------------------------------------------\n";
   */
  os.precision(oldPrec);

  return os;
}

struct GUTrack_v {
  int capacity;  // real number of tracks stored
  int numTracks; // real number of tracks stored
  int *status;   // status of the track: alive or killed (possible at rest ???)
  int *particleType;
  int *id;
  int *parentId; // index of the corresponding parent track in GeantTrack_v
  int *proc;     // index of physics process
  double *x;     // (x,y,z) position
  double *y;
  double *z;
  double *px; // momentum (px,py,pz)
  double *py;
  double *pz;
  double *E;      // total energy
  double *q;      // charge
  double *nint;   // number of interaction length left
  double *lambda; // interaction length
  double *s;      // step length
};

#endif
