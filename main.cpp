/*
 * 2次元FDTD法 Ver.003
 *
 * 説明：
 *  電子密度はIRIモデル、
 *  地磁気は IGRF モデル
 *
 * 必要なライブラリ：
 *  AndoLab, iri2016, igrf
 *
 * 出力先
 *  data_dir
 *
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>

#include <eigen3/Eigen/Dense>

#include <memory_allocate.h>
#include "fdtd2d.h"

std::string data_dir = "data/";

int main(int argc, char **argv){

  double ***Dr  = AndoLab::allocate_memory3d(2, Nr,   Nth+1, 0.0);
  double ***Dth = AndoLab::allocate_memory3d(2, Nr+1, Nth, 0.0);
  double ***Dph = AndoLab::allocate_memory3d(2, Nr+1, Nth+1, 0.0);
  double ***Er  = AndoLab::allocate_memory3d(2, Nr,   Nth+1, 0.0);;
  double ***Eth = AndoLab::allocate_memory3d(2, Nr+1, Nth, 0.0);
  double ***Eph = AndoLab::allocate_memory3d(2, Nr+1, Nth+1, 0.0);
  double **Hr   = AndoLab::allocate_memory2d(Nr+1, Nth, 0.0);
  double **Hth  = AndoLab::allocate_memory2d(Nr,   Nth+1, 0.0);
  double **Hph  = AndoLab::allocate_memory2d(Nr,   Nth, 0.0);
  
  /* Coefficient matrix to update E taking account into the ionosphere */
  Eigen::Matrix3d **C = new Eigen::Matrix3d* [Nr_iono];
  Eigen::Matrix3d **F = new Eigen::Matrix3d* [Nr_iono];
  Eigen::Matrix3d *C1 = new Eigen::Matrix3d [Nr_iono*(Nth+1)];
  Eigen::Matrix3d *F1 = new Eigen::Matrix3d [Nr_iono*(Nth+1)];

  for(int i = 0; i < Nr_iono; i++){
    C[i] = C1 + i*Nth;
    F[i] = F1 + i*Nth;
    for(int j = 0; j <= Nth; j++){
      C[i][j] << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
      F[i][j] << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    }
  }

  initialize_conductivity(C, F);

  /* PML only for +Theta-directed layer */
  double **Dr1    = AndoLab::allocate_memory2d(Nr,   PML_L+1, 0.0);
  double **Dr2    = AndoLab::allocate_memory2d(Nr,   PML_L+1, 0.0);
  double **Dph_r  = AndoLab::allocate_memory2d(Nr+1, PML_L+1, 0.0);
  double **Dph_th = AndoLab::allocate_memory2d(Nr+1, PML_L+1, 0.0);
  double **Hr1    = AndoLab::allocate_memory2d(Nr+1, PML_L,   0.0);
  double **Hr2    = AndoLab::allocate_memory2d(Nr+1, PML_L,   0.0);
  double **Hph_r  = AndoLab::allocate_memory2d(Nr,   PML_L,   0.0);
  double **Hph_th = AndoLab::allocate_memory2d(Nr,   PML_L,   0.0);
  
  double **Bph    = AndoLab::allocate_memory2d(2,   PML_L,   0.0);
  double *Bph_r   = new double [PML_L];
  double *Bph_th  = new double [PML_L];
  for(int i = 0; i < PML_L; i++){
    Bph_r[i] = 0.0;
    Bph_th[i] = 0.0;
  }

  double *C01 = new double [PML_L+1];
  double *C02 = new double [PML_L+1];
  double *C11 = new double [PML_L];
  double *C12 = new double [PML_L];
  initialize_pml(C01, C02, C11, C12);

  /* Vertical E-field at Earth's surface (f kHz) */
  std::complex <double> zj { 0., 1. };
  std::complex <double> *Er0 = new std::complex <double> [Nth+1 - PML_L];
  for(int i = 0; i <= Nth - PML_L; i++){
    Er0[i] = std::complex <double> {0., 0.};
  }

  /* Surface impedance */
  double *Rs = new double [Nth + 1];
  double *Ls = new double [Nth + 1];
  initialize_surface_impedance(Rs, Ls);

  ///時間ループ///
  for(int n = 1; n <= Nt; n++){
    if ( n%100 == 0 ){
      std::cout << n << " / " << Nt << "\n";
    }
    
    int NEW = (n+1) % 2;
    int OLD = n % 2;
    
    update_Dr(Dr, Hph, NEW, OLD);
    update_Dth(Dth, Hph, NEW, OLD);
    update_Dph(Dph, Hr, Hth, NEW, OLD);

    update_Dr_PML(Dr[NEW], Dr1, Dr2, Hph, C01, C02);
    update_Dph_PML(Dph[NEW], Dph_r, Dph_th, Hth, Hr, C01, C02);

    double t = (n - 0.5) * Dt;
    Dr[NEW][0][0] -= Dt * Jr(t); //θ=0, i=j=0  //Jr

    update_Er(Er, Eth, Eph, Dr, Dth, Dph, NEW, OLD, C, F);
    update_Eth(Er, Eth, Eph, Dr, Dth, Dph, NEW, OLD, C, F);
    update_Eph(Er, Eth, Eph, Dr, Dth, Dph, NEW, OLD, C, F);

    update_Hr(Hr, Eph[NEW]);
    update_Hth(Hth, Eph[NEW], Rs, Ls);
    update_Hph(Hph, Er[NEW], Eth[NEW], Rs, Ls);

    update_Hr_PML(Hr, Hr1, Hr2, Eph[NEW], C11, C12);
    update_Hph_PML(Hph, Hph_r, Hph_th, Er[NEW], Eth[NEW], C11, C12,
        Rs, Ls, Bph, Bph_r, Bph_th, NEW);

//    if ( n%10 == 0 ) output(Er, NEW, n);

    /* 地表面電界のフーリエ変換 */
    t = n * Dt;
    for(int j = 0; j <= Nth - PML_L; j++){
      Er0[j] += Er[NEW][0][j] * std::exp( -1.0 * zj * OMG * t ) * Dt;
    }
  }

  /* 地表面電界強度の出力 */
  std::ofstream ofs(data_dir + "Er_surface.dat");
  for(int j = 0; j <= Nth - PML_L; j++){
    ofs << j * R0 * dth * 1e-3 << " " << std::abs( Er0[j] ) << " "
        << std::arg( Er0[j] ) << "\n";
  }
  ofs.close();

  AndoLab::deallocate_memory3d(Dr);
  AndoLab::deallocate_memory3d(Dth);
  AndoLab::deallocate_memory3d(Dph);
  AndoLab::deallocate_memory3d(Er);
  AndoLab::deallocate_memory3d(Eth);
  AndoLab::deallocate_memory3d(Eph);
  AndoLab::deallocate_memory2d(Hr);
  AndoLab::deallocate_memory2d(Hth);
  AndoLab::deallocate_memory2d(Hph);
  delete [] C1;
  delete [] C;
  delete [] F1;
  delete [] F;
  delete [] Er0;

  AndoLab::deallocate_memory2d(Dr1);
  AndoLab::deallocate_memory2d(Dr2);
  AndoLab::deallocate_memory2d(Dph_r);
  AndoLab::deallocate_memory2d(Dph_th);
  AndoLab::deallocate_memory2d(Hr1);
  AndoLab::deallocate_memory2d(Hr2);
  AndoLab::deallocate_memory2d(Hph_r);
  AndoLab::deallocate_memory2d(Hph_th);
  AndoLab::deallocate_memory2d(Bph);

  delete [] Bph_r;
  delete [] Bph_th;
  delete [] C01;
  delete [] C02;
  delete [] C11;
  delete [] C12;

  delete [] Rs;
  delete [] Ls;

  return 0;
}
