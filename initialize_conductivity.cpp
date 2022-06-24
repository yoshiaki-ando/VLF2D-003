/*
 * initialize_conductivity.cpp
 *
 * 地理座標: 北極南極がz軸上、経度0°が +x軸上にあるものとする
 * 計算座標: 送信点を +z軸上、伝搬パスは zx平面内とする
 *
 *  Created on: 2021/06/02
 *      Author: ando
 */
#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <eigen3/Eigen/Dense>

#include <Vector3d.h>
#include <iri2016.h>
#include <geomag.h>

#include "fdtd2d.h"

extern std::string data_dir;

constexpr double M2KM { 1e-3 }; /* z[m] -> [km] */

///* Wait and Spies モデルによる電子密度 */
//double N(double z_in_m){
//  return 1.43e13 * exp( -0.15 * z_prime ) * exp( (beta-0.15) * (z_in_m*M2KM - z_prime) );
//}

/* 衝突周波数 */
double nu(const double z_in_m){
  return 4.303e11 * exp( -0.1622 * z_in_m*M2KM );
}

double omg_p(const double Ne){ /* Ne [m^-3] からプラズマ周波数を計算する */
  return sqrt( Ne*CHARGE_e*CHARGE_e/(MASS_e*EPS0) );
}

constexpr double PHI { 0.0 }; /* 2次元なので、シミュレーション座標でΦ = 0とする */
constexpr double s_ph { std::sin(PHI) };
constexpr double c_ph { std::cos(PHI) };

/*
 * 導電率テンソルの設定
 */
void initialize_conductivity(Eigen::Matrix3d** C, Eigen::Matrix3d** F){

  std::complex <double> zj{ 0.0, 1.0 }; /* 虚数単位 */

  /* 送受信点座標 */
  AndoLab::Vector3d <double> Tx =
      AndoLab::geographic_coordinate(Tx_Latitude, Tx_Longitude);
  AndoLab::Vector3d <double> Rx =
      AndoLab::geographic_coordinate(Rx_Latitude, Rx_Longitude);

  const double cosAlpha = Tx % Rx; /* 送受信点がなす角α */

  /* 計算座標の単位ベクトル */
  AndoLab::Vector3d <double> X = ( (Rx - cosAlpha * Tx) / std::sqrt(1 - cosAlpha*cosAlpha) ).n();
  AndoLab::Vector3d <double> Z =
        AndoLab::geographic_coordinate(Tx_Latitude, Tx_Longitude);
  AndoLab::Vector3d <double> Y = Z * X;

  /* 本来ここで設定するものではない */
  constexpr int Year { 2015 };
  constexpr int Month { 10 };
  constexpr int Day { 1 };

  AndoLab::iri2016 iri; /* IRI2016 */
  float *Ne = new float [Nr_iono];
  iri.set_datetime(Year, Month, Day, 17, 0); /* UT */

  std::string filename = data_dir + "B2.dat";
//  std::ofstream ofs(filename.c_str());

  std::ofstream ofs_ne( (data_dir + "Ne.dat").c_str() );
  for(int j = 0; j <= Nth; j++){
//  for(int j = 0; j <= Nth; j+=300){
    if (j%10 == 0) std::cout << "j = " << j << " / " << Nth << std::endl;

    const double THETA = j*dth;
    const double s_th = std::sin(THETA);
    const double c_th = std::cos(THETA);

    /* デカルト⇔球座標変換 */
    Eigen::Matrix3d Ph2Car, Car2Ph;
    Ph2Car <<
        s_th*c_ph , c_th*c_ph , -s_ph,
        s_th*s_ph , c_th*s_ph , c_ph ,
        c_th    , -s_th   , 0.0;
    Car2Ph <<
        s_th*c_ph , s_th*s_ph , c_th,
        c_th*c_ph , c_th*s_ph , -s_th,
        -s_ph   , c_ph    , 0.0;

    /* 現在の位置 */
    AndoLab::Vector3d <double> Pos = Z * std::cos(j*dth) + X * std::sin(j*dth);

    /* 現在の位置におけるθ方向 */
//    AndoLab::Vector3d <double> THETA = X * std::cos(j*dth) - Z * std::sin(j*dth);

    /* IRIを用いたオリジナル(IRI + Danilov)の電子密度 */
    iri.set_coord( float(Pos.latitude()), float(Pos.longitude()) );
    AndoLab::original_model(Nr_iono, float(Lower_boundary_of_ionosphere*M2KM),
        float(dr*M2KM), iri, Ne);

    for(int i = 0; i < Nr_iono; i++){ //Er(i+1/2,j)
      double z = Lower_boundary_of_ionosphere + i*dr;
      ofs_ne << j*dth*R0 * M2KM << " " << z*M2KM << " " << Ne[i] << std::endl;

      /* 地磁気 */
      double F0, Inc, Dec;

      /* IGRFモデルによる地磁気を取得 */
      calc_geomagnetic_field(
          Year, Month, Day,
          z*M2KM, Pos.latitude(), Pos.longitude(),
          Dec, Inc, F0);
      Dec *= AndoLab::DEG2RAD;
      Inc *= AndoLab::DEG2RAD;
      F0 *= 1e-9; /* [nT] to [T] */

      /* 地理座標での磁場 */
      AndoLab::Vector3d <double> B0_geo = F0 * (
          - std::sin(Inc) * Pos.r_vector()
          - std::cos(Inc)*std::cos(Dec) * Pos.theta_vector()
          + std::cos(Inc)*std::sin(Dec) * Pos.phi_vector() );

      /* 磁場の方向 */
      const double THE0 { std::acos(B0_geo.n()%Z) };
      const double PHI0 { std::atan2(B0_geo.n()%Y, B0_geo.n()%X) };
//      if (i == 0){
//        ofs << B0_geo.n()%X << " " << B0_geo.n()%Y << " " << B0_geo.n()%Z << " ";
//        ofs << 1.0 << " " << 0.0 << " " << 0.0 << " "
//            << 0.0 << " " << 1.0 << " " << 0.0 << " "
//            << 0.0 << " " << 0.0 << " " << 1.0 << " ";
//        ofs << B0_geo.n().x() << " " << B0_geo.n().y() << " " << B0_geo.n().z() << " ";
//        ofs << Pos.theta_vector().x() << " " << Pos.theta_vector().y() << " " << Pos.theta_vector().z() << " ";
//        ofs << Pos.phi_vector().x() << " " << Pos.phi_vector().y() << " " << Pos.phi_vector().z() << " ";
//        ofs << Dec * AndoLab::RAD2DEG << " " << Inc * AndoLab::RAD2DEG << " " << F0*1e9 << " ";
//        ofs << THE0 * AndoLab::RAD2DEG << " " << PHI0 * AndoLab::RAD2DEG << std::endl;
//      }
//      continue;

      /* 地磁気方向の変換 */
      Eigen::Matrix3d R1, R2;
      R1<< cos(THE0), 0.0, sin(THE0),
          0.0,        1.0, 0.0,
          -sin(THE0), 0.0, cos(THE0);
      R2<<cos(PHI0), -sin(PHI0), 0.0,
          sin(PHI0),  cos(PHI0), 0.0,
          0.0,        0.0,       1.0;

      double eps_r = Refractive_index(z) * Refractive_index(z);

      const double Omg_p { omg_p(Ne[i]) }; /* プラズマ周波数 */
      const double OMG_c { CHARGE_e * F0 / MASS_e }; /* サイクロトロン周波数 */

      std::complex <double> omg_prime { OMG - zj*nu(z) };
      std::complex <double> alpha { omg_prime / (OMG_c*OMG_c - omg_prime*omg_prime) };
      std::complex <double> beta { zj * OMG_c / (OMG_c*OMG_c - omg_prime*omg_prime) };
      Eigen::Matrix3cd SIGz;
      SIGz <<
          alpha, beta,  0.0,
          -beta, alpha, 0.0,
          0.0,   0.0,   -1.0/omg_prime;
      SIGz = (zj * EPS0 * Omg_p * Omg_p) * SIGz;

      Eigen::Matrix3d SIGcar = R2 * R1 * SIGz.real() * R1.inverse() * R2.inverse();
      Eigen::Matrix3d SIG = Car2Ph * SIGcar * Ph2Car;

      /////(A,B,)C,Fの設定/////
      Eigen::Matrix3d A;
      A = EPS0*eps_r/Dt * Eigen::Matrix3d::Identity() + 0.5 * SIG;
      Eigen::Matrix3d B;
      B = EPS0*eps_r/Dt * Eigen::Matrix3d::Identity() - 0.5 * SIG;

      C[i][j] = A.inverse() * B;
      F[i][j] = (1.0/Dt) * A.inverse();

    }

//    ofs << "\n";
    ofs_ne << "\n";
  }

//  ofs.close();
  ofs_ne.close();
//  exit(0);
}
