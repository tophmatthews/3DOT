#include "bubble.h"

bubbleBase::bubbleBase()
{
  _rad = simconf->bub_rad;
  if(simconf->bub_model > 20)
    _temp = 2000;
  else
    _temp = 1000;
  
  if ( simconf->bub_model == 0 )
    _model = CONSTANT;
  else if ( simconf->bub_model == 21 || simconf->bub_model == 11 )
    _model = VDW;
  else
    _model = RONCHI;
  
  calcArho();
}
  

void bubbleBase::calcArho( )
{
  double m = 131.3;
  switch( _model )
  {
    case CONSTANT:
      arho = 3.0; // atoms/m3
      break;
    case VDW:
    {
      double R = _rad / 1e10; // convert r from A to m units
      double B = 8.5E-29; // m^3/atom
      double k = 1.38E-23; // J/K: boltzmann's constant
      double gamma = 1; // N/A : surface tension
      double sigma = 0; // Pa : hydrostatic tension
      
      arho = 1 / ( B + 1 / ( ( 2.0*gamma/k/_temp )/R + sigma / k / _temp ) ); // atoms/m3
      break;
    }
    case RONCHI:
    { // Uses atomic density calculated from Ronchi JNM 96 1981 if above > 1000 Ang
      // Atomic density [atom/m3], and corresponding radius [ang]
      alglib::real_1d_array aden_1000 = "[55.450784356454,55.5377957334436,55.6331059132479,62.6370814939524,63.3302286745124,63.7356937826205,64.0239059372474,64.2462544706319,64.4288409631805,64.5832566489719,64.7159932342954,64.8351012999588,64.9407270324329,65.035241772715,65.1235792339727,65.2030912969004,65.2745502608825,65.3435431323695,65.4091404148553,65.4702948380086,65.5258646891634,65.5796673951564,65.6312252048213]";
      alglib::real_1d_array aden_2000 = "[54.757637175894,54.8446485528837,54.939958732688,62.6370814939524,63.3302286745124,63.7356937826205,64.0239059372474,64.2462544706319,64.4288409631805,64.5832566489719,64.7159932342954,64.8351012999588,64.9407270324329,65.035241772715,65.1235792339727,65.2030912969004,65.2745502608825,65.3435431323695,65.4091404148553,65.4702948380086,65.5258646891634,65.5796673951564,65.6312252048213]";
      alglib::real_1d_array rad_1000 = "[13.9978321147582,13.9108207377686,13.8155105579643,6.76799333660698,5.99547256850552,5.48224020470897,5.04911628091454,4.64339139880829,4.24540003714206,3.84881349203026,3.45428175275461,3.06433236504625,2.68165172724768,2.30774839995093,1.94302431707686,1.58691839100423,1.23775367362886,0.892867216930567,0.548415001603923,0.198914460767731,-0.163622819986412,-0.551416551592393,-0.984291377934894]";
      alglib::real_1d_array rad_2000 = "[13.9978321147582,13.9108207377686,13.8155105579643,6.02192375459269,5.21765946353058,4.68583808905555,4.25416126367599,3.87016155133323,3.5098967985855,3.16237735096294,2.82270376367163,2.48849307332342,2.15913431058908,1.83442691554746,1.51428683619517,1.19838842603204,0.885919153228741,0.575568610249621,0.265333668788019,-0.0476085139147253,-0.367254284844385,-0.699318101203609,-1.05214592238772]";
      
      alglib::spline1dinterpolant spline;
      if (_temp == 1000)
        alglib::spline1dbuildcubic(rad_1000, aden_1000, spline);
      else if (_temp == 2000)
        alglib::spline1dbuildcubic(rad_2000, aden_2000, spline);
      else
        fprintf( stderr, "\n!!!*** Unsupported temperature for the Ronchi Bubble EOS ***!!! \n Defaulting to 2000K\n");

      double lnr = log( _rad );
      arho = alglib::spline1dcalc(spline, lnr);
      arho = exp( arho );

      break;
    }
  }
  
  rho = arho / 6.022e23 * m / 1.0e6; // g/cc
  arho /= 1.0e30; // [atom/ang]

  std::cout << "bubble rho [g/cc]: " << rho << " arho [atoms/Ang^3]: " << arho << endl;
}