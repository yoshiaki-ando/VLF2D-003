#include "fdtd2d.h"

/* 電流モーメント I Δl */
double I_Dl(double t){
  return (t-t0) / s * exp( - (t-t0)*(t-t0) / (2.0*s*s) );
}

double Jr(double t){   ///Jr(R0+(i+1/2)*Dr,j*Dth) = Jr(i+1/2,j)
  double R = (R0 + dr/2.0) * dth/2.0; /* radius of circular area at current density, 地表面に固定 */
  return I_Dl(t) / dr / M_PI / R / R; /* I Δl / ΔV */
}
