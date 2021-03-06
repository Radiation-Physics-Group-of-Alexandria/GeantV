/*
 * ConstFieldHelixStepper.h
 *
 *  Created on: Apr 23, 2014
 *      Author: swenzel
 */

#ifndef CONSTFIELDHELIXSTEPPER_H_
#define CONSTFIELDHELIXSTEPPER_H_

#include "Geant/Config.h"

// add the sincos function on MAC because sincos is not part of math.h
#if defined(__APPLE__) & !defined(NVCC)  // possibly other conditions
inline
void sincos(double x, double *s, double *c) {
  __sincos(x,s,c);
}
#endif

namespace Geant
{

  /**
  * A very simple stepper treating the propagation of particles in a constant Bz magnetic field
  * ( neglecting energy loss of particle )
  * This class is roughly equivalent to TGeoHelix in ROOT
  */
  class ConstBzFieldHelixStepper
  {
    private:
      double fBz;

    public:
      VECCORE_ATT_HOST_DEVICE
      ConstBzFieldHelixStepper( double Bz = 0. ) : fBz(Bz) {}

      void SetBz( double Bz ){ fBz = Bz; }
      double GetBz() const { return fBz; }

      /*
      template<typename RT, typename Vector3D>
      RT GetCurvature( Vector3D const & dir,
              double const charge, double const momentum ) const
      {
        if (charge == 0) return RT(0.);
        return abs( kB2C * fBz * dir.FastInverseScaledXYLength( momentum ) );
      }
*/

      /**
       * this function propagates the track along the helix solution by a step
       * input: current position, current direction, some particle properties
       * output: new position, new direction of particle
       */
       template<typename BaseType, typename BaseIType>
       inline
       __attribute__((always_inline))
       VECCORE_ATT_HOST_DEVICE
       void DoStep( BaseType const & /*posx*/, BaseType const & /*posy*/, BaseType const & /*posz*/,
                    BaseType const & /*dirx*/, BaseType const & /*diry*/, BaseType const & /*dirz*/,
                    BaseIType const & /*charge*/, BaseType const & /*momentum*/, BaseType const & /*step*/,
                    BaseType & /*newsposx*/, BaseType  & /*newposy*/, BaseType  & /*newposz*/,
                    BaseType & /*newdirx*/, BaseType  & /*newdiry*/, BaseType  & /*newdirz*/
                  ) const ;

       /**
        * basket version of dostep
        * version that takes plain arrays as input; suited for current Geant-V
	* 
	* SW: for the moment (12.5.2015) commenting this out as not used 
	*
	*/
       //void DoStep_v( double const * /*posx*/, double const * /*posy*/, double const * /*posz*/,
       //               double const * /*dirx*/, double const * /*diry*/, double const * /*dirz*/,
       //              int const * /*charge*/, double const * /*momentum*/, double const * /*step*/,
       //              double * /*newsposx*/, double * /*newposy*/, double * /*newposz*/,
       //              double * /*newdirx*/, double * /*newdiry*/, double * /*newdirz*/,
       //              int np
       //           ) const ;

        // in future will offer versions that take containers as input

        /**
         * this function propagates the track along the helix solution by a step
         * input: current position, current direction, some particle properties
         * output: new position, new direction of particle
         */
         template<typename Vector3D, typename BaseType, typename BaseIType>
         VECCORE_ATT_HOST_DEVICE
         void DoStep( Vector3D const & pos, Vector3D const & dir,
                           BaseIType const & charge, BaseType const &  momentum,
                           BaseType const & step,
                           Vector3D & newpos,
                           Vector3D & newdir
                         ) const
         {
             DoStep(pos[0],pos[1],pos[2],dir[0],dir[1],dir[2], charge, momentum, step,
                         newpos[0], newpos[1], newpos[2], newdir[0], newdir[1], newdir[2]);
         }

  }; // end class declaration

  /**
   * this function propagates the track along the "helix-solution" by a step step
   * input: current position (x0, y0, z0), current direction ( dx0, dy0, dz0 ), some particle properties
   * output: new position, new direction of particle
   */
  template<typename BaseDType, typename BaseIType>
  inline
  __attribute__((always_inline))
  void ConstBzFieldHelixStepper::DoStep(
               BaseDType const & x0, BaseDType const & y0, BaseDType const & z0,
               BaseDType const & dx0, BaseDType const & dy0, BaseDType const & dz0,
               BaseIType const & charge, BaseDType const & momentum, BaseDType const & step,
               BaseDType & x, BaseDType & y, BaseDType & z,
               BaseDType & dx, BaseDType & dy, BaseDType & dz
             ) const
  {
      const double kB2C_local = -0.299792458e-3;
      const double kSmall = 1.E-30;
      // could do a fast square root here
      BaseDType dt = sqrt((dx0*dx0) + (dy0*dy0)) + kSmall;
      BaseDType invnorm=1./dt;
      // radius has sign and determines the sense of rotation
      BaseDType R = momentum*dt/((kB2C_local*BaseDType(charge))*(fBz));

      BaseDType cosa= dx0*invnorm;
      BaseDType sina= dy0*invnorm;
      BaseDType phi = step * BaseDType(charge) * fBz * kB2C_local / momentum;

      BaseDType cosphi;
      BaseDType sinphi;
      sincos(phi, &sinphi, &cosphi);

      x = x0 + R*( -sina - ( -cosphi*sina - sinphi*cosa ));
      y = y0 + R*( cosa  - ( -sinphi*sina + cosphi*cosa ));
      z = z0 + step * dz0;

      dx = dx0 * cosphi - sinphi * dy0;
      dy = dx0 * sinphi + cosphi * dy0;
      dz = dz0;
  }

  /**
   * basket version of dostep
   */
  /*
   SW: commented out due to explicit Vc dependence and since it is not currently used
       leaving the code here to show how one would dispatch to the kernel with Vc
#define _R_ __restrict__
  void ConstBzFieldHelixStepper::DoStep_v(
                        double const * _R_ posx, double const * _R_ posy, double const * _R_ posz,
                        double const * _R_ dirx, double const * _R_ diry, double const * _R_ dirz,
                        int const * _R_ charge, double const * _R_ momentum, double const * _R_ step,
                        double * _R_ newposx, double * _R_ newposy, double * _R_ newposz,
                        double * _R_ newdirx, double * _R_ newdiry, double * _R_ newdirz,
                        int np
                     ) const
   {
       // we have choice here: ( try autovectorization: )

//#pragma ivdep
//      for (int i=0;i<np;++i){
//            DoStep( posx[i], posy[i], posz[i], dirx[i], diry[i], dirz[i],
//                    charge[i], momentum[i], step[i],
//                    newposx[i], newposy[i], newposz[i],
//                    newdirx[i], newdiry[i], newdirz[i]
//                  );
//       }

       // alternative loop with Vc:

       typedef Vc::Vector<double> Vcdouble_t;
       typedef Vc::Vector<int> Vcint_t;
       for (int i=0;i<np;i+=Vcdouble_t::Size)
       {
            // results cannot not be temporaries
            Vcdouble_t newposx_v, newposy_v, newposz_v,
                       newdirx_v, newdiry_v,newdirz_v;
            DoStep( Vcdouble_t(posx[i]), Vcdouble_t(posy[i]), Vcdouble_t(posz[i]),
                    Vcdouble_t(dirx[i]), Vcdouble_t(diry[i]), Vcdouble_t(dirz[i]),
                    Vcint_t(charge[i]), Vcdouble_t(momentum[i]), Vcdouble_t(step[i]),
                    newposx_v,
                    newposy_v,
                    newposz_v,
                    newdirx_v,
                    newdiry_v,
                    newdirz_v
                   );
            // write results
            newposx_v.store(&newposx[i]);
            newposy_v.store(&newposy[i]);
            newposz_v.store(&newposz[i]);
            newdirx_v.store(&newdirx[i]);
            newdiry_v.store(&newdiry[i]);
            newdirz_v.store(&newdirz[i]);
       }
       // tail part: tobedone
   }
  */

   //TODO: above stepper is tailored/specialized to B=(0,0,Bz) in the global frame of reference
   // might need to provide more general class in which the constant field has arbitrary direction


} // end geant namespace

#endif /* CONSTFIELDHELIXSTEPPER_H_ */
