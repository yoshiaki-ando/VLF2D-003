#include <cmath>
#include <eigen3/Eigen/Dense>

/* Physical constants */
constexpr double C0 { 3.0e8 };
constexpr double MU0 { 4.0 * M_PI * 1e-7 };
constexpr double EPS0 { 1.0 / MU0 / C0 / C0 };
constexpr double Z0 { std::sqrt( MU0 / EPS0 ) };

constexpr double CHARGE_e { 1.602e-19 }; //[C]
constexpr double MASS_e { 9.1e-31 }; //[kg]

constexpr double R0 { 6370.0e3 }; /* Radius of Earth */

/************************************************************
 * Setting parameters (from here)
 ************************************************************/
constexpr double Rr  { 120.0e3 }; /* Analysis Region in r-direction */
constexpr double Rth { 2200.0e3 }; /* Analysis Region in θ-direction */
constexpr double dr { 0.5e3 };  /* Cell size in r */
constexpr double Rdth { 0.5e3 }; /* Cell size in θ */
constexpr double Tmax { 0.00848528 };

constexpr double FREQ { 22.2e3 }; /* [Hz] */
//constexpr double FREQ { 40.0e3 }; /* [Hz] */

/* えびの 緯度経度 */
constexpr double Tx_Latitude { 32.067650 };
constexpr double Tx_Longitude { 130.828978 };

/* おおたかどや山 緯度経度 */
//constexpr double Tx_Latitude { 37.372097972301226 };
//constexpr double Tx_Longitude { 140.84913351739013 };

/* 調布 緯度経度 */
constexpr double Rx_Latitude { 35.65966 };
constexpr double Rx_Longitude { 139.54328 };

/* 新潟 緯度経度 */
//constexpr double Rx_Latitude { 37.9 };
//constexpr double Rx_Longitude { 139.0 };

///* Daytime condition */
//constexpr double z_prime { 73.0 }; //[km]
//constexpr double beta { 0.3 };

///* Nighttime condition */
//constexpr double z_prime { 85.0 }; //[km]
//constexpr double beta { 0.63 };

/* Current source Pulse waveform */
constexpr double s  { 1.0/2.0/M_PI/FREQ };
constexpr double t0 { 6.0 * s   };

/* PML */
constexpr int PML_L { 10 };
constexpr double PML_M { 3.2 };
constexpr double Gamma { 1.0e-6 };

/* Parameters of ionosphere */
/* この高度から電離圏の導電率を考慮する */
constexpr double Lower_boundary_of_ionosphere { 60.0e3 };

/************************************************************
 * Setting parameters (up to here)
 ************************************************************/


constexpr double dth { Rdth/R0 }; /* Angle making a cell */
constexpr int Nr { int(Rr/dr) };  /* Cell number in r */
constexpr int Nth { int(Rth/Rdth) }; /* Cell number in theta */

constexpr double Dt { 0.9 / C0 / sqrt( 1/(dr*dr) + 1/((R0*dth)*(R0*dth)) ) };
constexpr int Nt { int(Tmax/Dt)+1 };

constexpr double OMG { FREQ * 2. * M_PI };

/* Parameters of ionosphere */
constexpr int Nr_iono { int( (Rr - Lower_boundary_of_ionosphere) / dr ) + 1 };
constexpr int Nr_atmo { int( Lower_boundary_of_ionosphere / dr ) };

double I_Dl(double t);
double Jr(double t);

Eigen::Matrix3d define_SIGzr(double z);

void update_Dr(double ***Dr, double **Hph, const int NEW, const int OLD);
void update_Dth(double ***Dth, double **Hph, const int NEW, const int OLD);
void update_Dph(double ***Dph, double **Hr, double **Hth, const int NEW, const int OLD);
void update_Dr_PML(double **Dr, double **Dr1, double **Dr2,
    double **Hph, double *C1, double *C2);
void update_Dph_PML(double **Dph, double **Dph_r, double **Dph_th,
    double **Hth, double **Hr,
    double *C1, double *C2);

void update_Er(double ***Er, double ***Eth, double ***Eph,
    double ***Dr, double ***Dth, double ***Dph,
    const int NEW, const int OLD,
    Eigen::Matrix3d **C, Eigen::Matrix3d **F);
void update_Eth(double ***Er, double ***Eth, double ***Eph,
    double ***Dr, double ***Dth, double ***Dph,
    const int NEW, const int OLD,
    Eigen::Matrix3d **C, Eigen::Matrix3d **F);
void update_Eph(double ***Er, double ***Eth, double ***Eph,
    double ***Dr, double ***Dth, double ***Dph,
    const int NEW, const int OLD,
    Eigen::Matrix3d **C, Eigen::Matrix3d **F);

void update_Hr(double **Hr, double **Eph);
void update_Hth(double **Hth, double **Eph, double *Rs, double *Ls);
void update_Hph(double **Hph, double **Er, double **Eth, double *Rs, double *Ls);
void update_Hr_PML(double **Hr, double **Hr1, double **Hr2,
    double **Eph, double *C1, double *C2);
void update_Hph_PML(double **Hph, double **Hph_r, double **Hph_th,
    double **Er, double **Eth, double *C1, double *C2,
    double *Rs, double *Ls, double **Bph, double *Bph_r, double *Bph_th, const int NEW);

void initialize_conductivity(Eigen::Matrix3d** C, Eigen::Matrix3d** F);

void initialize_pml(double *C01, double *C02, double *C11, double *C12);

void initialize_surface_impedance(double *Rs, double *Ls);

void set_perturbation_parameter(int argc, char **argv,
    double &Lp, double &z_dec, double &sig_per);

std::string suffix(double Lp, double z_dec, double sig_per);

void output(double ***Er, int NEW, int n);

inline double r(double i){
  return R0 + i*dr;
}

inline double Refractive_index(const double z){
  /* refractive index of standard air */
  return 1.000325 - 0.039e-6 * z;
}
