#include "bubble.h"

bubbleBase::bubbleBase()
{
  _rad = simconf->bub_rad;
  if(simconf->bub_model > 20)
    _temp = 2000;
  else
    _temp = 1000;
  
  if ( simconf->bub_model == 21 || simconf->bub_model == 11 || _rad > 1000 )
    _model = VDW;
  else
    _model = RONCHI;
  
  calcArho();
}
  

void bubbleBase::calcArho( )
{
  double R = _rad / 1e10; // convert r from A to m units
  double B = 8.5E-29; // m^3/atom
  double k = 1.38E-23; // J/K: boltzmann's constant
  double gamma = 1; // N/A : surface tension
  double sigma = 0; // Pa : hydrostatic tension
  double m = 132.0;
  switch( _model )
  {
    case VDW:
    {
      arho = 1 / ( B + 1 / ( ( 2*gamma/k/_temp )/R + sigma / k / _temp ) ); // atoms/m3
      break;
    }
    case RONCHI:
    { // Uses atomic density calculated from Ronchi JNM 96 1981 if above > 1000 Ang
      // Atomic density [atom/m3], and corresponding radius [ang]
      alglib::real_1d_array aden_1000 = "[1.29E+27,1.43E+26,1.45E+25,1.45E+24,1.60E+27,3.19E+27,4.79E+27,6.39E+27,7.98E+27,9.57E+27,1.12E+28,1.28E+28,1.44E+28,1.60E+28,1.76E+28,1.92E+28,2.08E+28,2.23E+28,2.39E+28,2.55E+28,2.71E+28,2.87E+28,3.03E+28,3.19E+28]";
      alglib::real_1d_array aden_2000 = "[6.83E+26,7.20E+25,7.24E+24,7.25E+23,1.60E+27,3.19E+27,4.79E+27,6.39E+27,7.98E+27,9.57E+27,1.12E+28,1.28E+28,1.44E+28,1.60E+28,1.76E+28,1.92E+28,2.08E+28,2.23E+28,2.39E+28,2.55E+28,2.71E+28,2.87E+28,3.03E+28,3.19E+28]";
      alglib::real_1d_array rad_1000 = "[1000,10000,100000,1000000,869.5652174,401.6064257,240.3846154,155.8846454,103.8961039,69.78367062,46.93733865,31.63555837,21.42015637,14.6092038,10.0517666,6.979828296,4.888660751,3.447859741,2.442121715,1.730507991,1.220077597,0.849062211,0.57613311,0.373703948]";
      alglib::real_1d_array rad_2000 = "[1000,10000,100000,1000000,412.371134,184.501845,108.401084,70.39774727,47.95013186,33.44481605,23.62669817,16.82227269,12.04311435,8.663634395,6.261544723,4.546177801,3.314770618,2.425212509,1.778141309,1.303865963,0.953506999,0.692633496,0.49692404,0.349187615]";

      alglib::spline1dinterpolant spline;
      if (_temp == 1000)
        alglib::spline1dbuildcubic(rad_1000, aden_1000, spline);
      else if (_temp == 2000)
        alglib::spline1dbuildcubic(rad_2000, aden_2000, spline);
      else
        fprintf( stderr, "\n!!!*** Unsupported temperature for the Ronchi Bubble EOS ***!!! \n Defaulting to 2000K\n");

      arho = alglib::spline1dcalc(spline, _rad);

      break;
    }
  }
  rho = arho / 6.022e23 * m / 1e6; // g/cc
  arho /= 1e30; // [atom/ang]

  std::cout << "rho: " << rho << " arho: " << arho << endl;
}