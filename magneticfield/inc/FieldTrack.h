#ifndef _FIELDTRACK_H_
#define _FIELDTRACK_H_

struct FieldTrack{

public:
  //data members

  double PosMomVector[6];

  //And functions 
  void DumpToArray(double valArr[]){ //12 from ncompSVEC as in both TemplateGUIntegrationDriver
    for (int i = 0; i < 6; ++i)        //and GUFieldTrack function
    {
      valArr[i] = PosMomVector[i];
    }
  }

  void LoadFromArray(const double valArr[], int noVarsIntegrated = 6){
    for (int i = 0; i < noVarsIntegrated; ++i)
    {
      PosMomVector[i] = valArr[i];
    }
  }

  void SetCurveLength(double len){
    fDistanceAlongCurve = len;
  }

  double GetCurveLength(){
    return fDistanceAlongCurve;
  }

private: 
  double fDistanceAlongCurve = 0.0;

};



#endif 