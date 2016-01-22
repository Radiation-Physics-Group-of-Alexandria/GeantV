

#ifndef TVecUniMF_H
#define TVecUniMF_H

#include <iostream>
#include "base/Vector3D.h"


class TVecUniMF
{  
    public:  

        TVecUniMF()
        {
          std::cout<<"-- entered TVecUniMF constructor ---"<<std::endl;
        }

        ~TVecUniMF() {}

    private:
        vecgeom::Vector3D<typename Vc::Vector<float>> fFieldComponents;

};

#endif
