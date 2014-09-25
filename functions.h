#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <vector>

inline void v_cross( const double *a1, const double *a2, double *b )
{
    for( int i = 0; i < 3; i++ )
      b[i] = a1[(i+1)%3] * a2[(i+2)%3] - a1[(i+2)%3] * a2[(i+1)%3];
}

inline void v_scale( double *a1, double b )
{
    for( int i = 0; i < 3; i++ ) a1[i] = a1[i] * b;
}

inline double v_dot( const double *a1, const double *a2 )
{
    double b = 0.0;
      for( int i = 0; i < 3; i++ ) b += a1[i] * a2[i];
        return b;
}

inline void v_norm( double *a1, double b = 1.0 ) // in-place normalize to b ( = 1.0 default )
{
    v_scale( a1, b / sqrt( v_dot( a1, a1) ) );
}

inline double sqr( double a ) { return a*a; }
inline double cub( double a ) { return a*a*a; }

// random numbers
const double drm = double(RAND_MAX)+1.0;
inline double dr250() { return double(rand())/drm; };
inline void r250_init( int s ) { srand(s); }

// Calculate atom density in bubble
inline double calc_rho( double r, double m )
{
  double R = r / 1.0e10; // convert r from A to m units
  double B = 8.5E-29; // m^3/atom
  double k = 1.38E-23; // J/K: boltzmann's constant
  double T = 1000; // K : temperature
  double gamma = 1; // N/A : surface tension
  double sigma = 0; // Pa : hydrostatic tension
  double arho; // atomic density
  double rho; // mass density
  
  arho = 1.0 / ( B + 1.0 / ( ( 2.0*gamma/k/T )/R + sigma / k / T ) ); // atoms/m3
  rho = arho / 6.022e23 * m / 1e6; // g/cc
  //std::cout << "rho: " << rho << " arho: " << arho << endl;
  
  return rho;
}

inline bool inbubble( double bub_cen[3], double pos[3], double r )
{
  double a;
  for( int i = 0; i < 3; i++ )
    a += sqr(bub_cen[i]/2.0 - pos[i]);
  a = sqrt(a);
  if( a < r ) return true;
  else return false;
}

inline double avg( std::vector<int> x )
{
  double a = 0;
  for( int i=0; i < x.size(); i++ )
    a += x.at(i);
  return a / x.size();
}

inline double avg( std::vector<long int> x )
{
  double a = 0;
  for( int i=0; i < x.size(); i++ )
    a += x.at(i);
  return a / x.size();
}

inline double stdev( std::vector<int> x )
{
  double a = avg(x);
  double b = 0;
  for( int i=0; i < x.size(); i++)
    b += sqr( x.at(i) - a );
  return sqrt( b / (x.size() - 1) );
}

inline double stdev( std::vector<long int> x )
{
  double a = avg(x);
  double b = 0;
  for( int i=0; i < x.size(); i++)
    b += sqr( x.at(i) - a );
  return sqrt( b / (x.size() - 1) );
}

inline double v_avg( std::vector<double> v)
{
  double v_avg = 0;
  for(std::vector<double>::iterator j=v.begin();j!=v.end();++j) v_avg += *j;
  v_avg /= v.size();
  return v_avg;
}

#endif

