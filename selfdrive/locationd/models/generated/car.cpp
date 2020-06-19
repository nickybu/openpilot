
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8228335980912175971) {
   out_8228335980912175971[0] = delta_x[0] + nom_x[0];
   out_8228335980912175971[1] = delta_x[1] + nom_x[1];
   out_8228335980912175971[2] = delta_x[2] + nom_x[2];
   out_8228335980912175971[3] = delta_x[3] + nom_x[3];
   out_8228335980912175971[4] = delta_x[4] + nom_x[4];
   out_8228335980912175971[5] = delta_x[5] + nom_x[5];
   out_8228335980912175971[6] = delta_x[6] + nom_x[6];
   out_8228335980912175971[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_7811466768422928389) {
   out_7811466768422928389[0] = -nom_x[0] + true_x[0];
   out_7811466768422928389[1] = -nom_x[1] + true_x[1];
   out_7811466768422928389[2] = -nom_x[2] + true_x[2];
   out_7811466768422928389[3] = -nom_x[3] + true_x[3];
   out_7811466768422928389[4] = -nom_x[4] + true_x[4];
   out_7811466768422928389[5] = -nom_x[5] + true_x[5];
   out_7811466768422928389[6] = -nom_x[6] + true_x[6];
   out_7811466768422928389[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_4345067464258930640) {
   out_4345067464258930640[0] = 1.0;
   out_4345067464258930640[1] = 0.0;
   out_4345067464258930640[2] = 0.0;
   out_4345067464258930640[3] = 0.0;
   out_4345067464258930640[4] = 0.0;
   out_4345067464258930640[5] = 0.0;
   out_4345067464258930640[6] = 0.0;
   out_4345067464258930640[7] = 0.0;
   out_4345067464258930640[8] = 0.0;
   out_4345067464258930640[9] = 1.0;
   out_4345067464258930640[10] = 0.0;
   out_4345067464258930640[11] = 0.0;
   out_4345067464258930640[12] = 0.0;
   out_4345067464258930640[13] = 0.0;
   out_4345067464258930640[14] = 0.0;
   out_4345067464258930640[15] = 0.0;
   out_4345067464258930640[16] = 0.0;
   out_4345067464258930640[17] = 0.0;
   out_4345067464258930640[18] = 1.0;
   out_4345067464258930640[19] = 0.0;
   out_4345067464258930640[20] = 0.0;
   out_4345067464258930640[21] = 0.0;
   out_4345067464258930640[22] = 0.0;
   out_4345067464258930640[23] = 0.0;
   out_4345067464258930640[24] = 0.0;
   out_4345067464258930640[25] = 0.0;
   out_4345067464258930640[26] = 0.0;
   out_4345067464258930640[27] = 1.0;
   out_4345067464258930640[28] = 0.0;
   out_4345067464258930640[29] = 0.0;
   out_4345067464258930640[30] = 0.0;
   out_4345067464258930640[31] = 0.0;
   out_4345067464258930640[32] = 0.0;
   out_4345067464258930640[33] = 0.0;
   out_4345067464258930640[34] = 0.0;
   out_4345067464258930640[35] = 0.0;
   out_4345067464258930640[36] = 1.0;
   out_4345067464258930640[37] = 0.0;
   out_4345067464258930640[38] = 0.0;
   out_4345067464258930640[39] = 0.0;
   out_4345067464258930640[40] = 0.0;
   out_4345067464258930640[41] = 0.0;
   out_4345067464258930640[42] = 0.0;
   out_4345067464258930640[43] = 0.0;
   out_4345067464258930640[44] = 0.0;
   out_4345067464258930640[45] = 1.0;
   out_4345067464258930640[46] = 0.0;
   out_4345067464258930640[47] = 0.0;
   out_4345067464258930640[48] = 0.0;
   out_4345067464258930640[49] = 0.0;
   out_4345067464258930640[50] = 0.0;
   out_4345067464258930640[51] = 0.0;
   out_4345067464258930640[52] = 0.0;
   out_4345067464258930640[53] = 0.0;
   out_4345067464258930640[54] = 1.0;
   out_4345067464258930640[55] = 0.0;
   out_4345067464258930640[56] = 0.0;
   out_4345067464258930640[57] = 0.0;
   out_4345067464258930640[58] = 0.0;
   out_4345067464258930640[59] = 0.0;
   out_4345067464258930640[60] = 0.0;
   out_4345067464258930640[61] = 0.0;
   out_4345067464258930640[62] = 0.0;
   out_4345067464258930640[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_5839489982569705744) {
   out_5839489982569705744[0] = state[0];
   out_5839489982569705744[1] = state[1];
   out_5839489982569705744[2] = state[2];
   out_5839489982569705744[3] = state[3];
   out_5839489982569705744[4] = state[4];
   out_5839489982569705744[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5839489982569705744[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5839489982569705744[7] = state[7];
}
void F_fun(double *state, double dt, double *out_5476922234156819524) {
   out_5476922234156819524[0] = 1;
   out_5476922234156819524[1] = 0;
   out_5476922234156819524[2] = 0;
   out_5476922234156819524[3] = 0;
   out_5476922234156819524[4] = 0;
   out_5476922234156819524[5] = 0;
   out_5476922234156819524[6] = 0;
   out_5476922234156819524[7] = 0;
   out_5476922234156819524[8] = 0;
   out_5476922234156819524[9] = 1;
   out_5476922234156819524[10] = 0;
   out_5476922234156819524[11] = 0;
   out_5476922234156819524[12] = 0;
   out_5476922234156819524[13] = 0;
   out_5476922234156819524[14] = 0;
   out_5476922234156819524[15] = 0;
   out_5476922234156819524[16] = 0;
   out_5476922234156819524[17] = 0;
   out_5476922234156819524[18] = 1;
   out_5476922234156819524[19] = 0;
   out_5476922234156819524[20] = 0;
   out_5476922234156819524[21] = 0;
   out_5476922234156819524[22] = 0;
   out_5476922234156819524[23] = 0;
   out_5476922234156819524[24] = 0;
   out_5476922234156819524[25] = 0;
   out_5476922234156819524[26] = 0;
   out_5476922234156819524[27] = 1;
   out_5476922234156819524[28] = 0;
   out_5476922234156819524[29] = 0;
   out_5476922234156819524[30] = 0;
   out_5476922234156819524[31] = 0;
   out_5476922234156819524[32] = 0;
   out_5476922234156819524[33] = 0;
   out_5476922234156819524[34] = 0;
   out_5476922234156819524[35] = 0;
   out_5476922234156819524[36] = 1;
   out_5476922234156819524[37] = 0;
   out_5476922234156819524[38] = 0;
   out_5476922234156819524[39] = 0;
   out_5476922234156819524[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_5476922234156819524[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_5476922234156819524[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5476922234156819524[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5476922234156819524[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_5476922234156819524[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_5476922234156819524[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_5476922234156819524[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_5476922234156819524[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_5476922234156819524[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_5476922234156819524[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5476922234156819524[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5476922234156819524[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_5476922234156819524[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_5476922234156819524[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_5476922234156819524[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5476922234156819524[56] = 0;
   out_5476922234156819524[57] = 0;
   out_5476922234156819524[58] = 0;
   out_5476922234156819524[59] = 0;
   out_5476922234156819524[60] = 0;
   out_5476922234156819524[61] = 0;
   out_5476922234156819524[62] = 0;
   out_5476922234156819524[63] = 1;
}
void h_25(double *state, double *unused, double *out_1729658063121940243) {
   out_1729658063121940243[0] = state[6];
}
void H_25(double *state, double *unused, double *out_5551778729974662920) {
   out_5551778729974662920[0] = 0;
   out_5551778729974662920[1] = 0;
   out_5551778729974662920[2] = 0;
   out_5551778729974662920[3] = 0;
   out_5551778729974662920[4] = 0;
   out_5551778729974662920[5] = 0;
   out_5551778729974662920[6] = 1;
   out_5551778729974662920[7] = 0;
}
void h_24(double *state, double *unused, double *out_1835918340960958870) {
   out_1835918340960958870[0] = state[4];
   out_1835918340960958870[1] = state[5];
}
void H_24(double *state, double *unused, double *out_5256232666140554348) {
   out_5256232666140554348[0] = 0;
   out_5256232666140554348[1] = 0;
   out_5256232666140554348[2] = 0;
   out_5256232666140554348[3] = 0;
   out_5256232666140554348[4] = 1;
   out_5256232666140554348[5] = 0;
   out_5256232666140554348[6] = 0;
   out_5256232666140554348[7] = 0;
   out_5256232666140554348[8] = 0;
   out_5256232666140554348[9] = 0;
   out_5256232666140554348[10] = 0;
   out_5256232666140554348[11] = 0;
   out_5256232666140554348[12] = 0;
   out_5256232666140554348[13] = 1;
   out_5256232666140554348[14] = 0;
   out_5256232666140554348[15] = 0;
}
void h_30(double *state, double *unused, double *out_5477911081103537199) {
   out_5477911081103537199[0] = state[4];
}
void H_30(double *state, double *unused, double *out_280600472658078466) {
   out_280600472658078466[0] = 0;
   out_280600472658078466[1] = 0;
   out_280600472658078466[2] = 0;
   out_280600472658078466[3] = 0;
   out_280600472658078466[4] = 1;
   out_280600472658078466[5] = 0;
   out_280600472658078466[6] = 0;
   out_280600472658078466[7] = 0;
}
void h_26(double *state, double *unused, double *out_5382044204511300316) {
   out_5382044204511300316[0] = state[7];
}
void H_26(double *state, double *unused, double *out_1277297089077403520) {
   out_1277297089077403520[0] = 0;
   out_1277297089077403520[1] = 0;
   out_1277297089077403520[2] = 0;
   out_1277297089077403520[3] = 0;
   out_1277297089077403520[4] = 0;
   out_1277297089077403520[5] = 0;
   out_1277297089077403520[6] = 0;
   out_1277297089077403520[7] = 1;
}
void h_27(double *state, double *unused, double *out_4755968270779530176) {
   out_4755968270779530176[0] = state[3];
}
void H_27(double *state, double *unused, double *out_5790801706628167388) {
   out_5790801706628167388[0] = 0;
   out_5790801706628167388[1] = 0;
   out_5790801706628167388[2] = 0;
   out_5790801706628167388[3] = 1;
   out_5790801706628167388[4] = 0;
   out_5790801706628167388[5] = 0;
   out_5790801706628167388[6] = 0;
   out_5790801706628167388[7] = 0;
}
void h_29(double *state, double *unused, double *out_2714519565818744062) {
   out_2714519565818744062[0] = state[1];
}
void H_29(double *state, double *unused, double *out_5984792589090309848) {
   out_5984792589090309848[0] = 0;
   out_5984792589090309848[1] = 1;
   out_5984792589090309848[2] = 0;
   out_5984792589090309848[3] = 0;
   out_5984792589090309848[4] = 0;
   out_5984792589090309848[5] = 0;
   out_5984792589090309848[6] = 0;
   out_5984792589090309848[7] = 0;
}
void h_28(double *state, double *unused, double *out_3466311643217169430) {
   out_3466311643217169430[0] = state[5];
   out_3466311643217169430[1] = state[6];
}
void H_28(double *state, double *unused, double *out_8659147327246331634) {
   out_8659147327246331634[0] = 0;
   out_8659147327246331634[1] = 0;
   out_8659147327246331634[2] = 0;
   out_8659147327246331634[3] = 0;
   out_8659147327246331634[4] = 0;
   out_8659147327246331634[5] = 1;
   out_8659147327246331634[6] = 0;
   out_8659147327246331634[7] = 0;
   out_8659147327246331634[8] = 0;
   out_8659147327246331634[9] = 0;
   out_8659147327246331634[10] = 0;
   out_8659147327246331634[11] = 0;
   out_8659147327246331634[12] = 0;
   out_8659147327246331634[13] = 0;
   out_8659147327246331634[14] = 1;
   out_8659147327246331634[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;
  
  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);
  
  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H); 
  
  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();
   

    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;
  
  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);
 
  // update cov 
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
