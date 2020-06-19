/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8228335980912175971);
void inv_err_fun(double *nom_x, double *true_x, double *out_7811466768422928389);
void H_mod_fun(double *state, double *out_4345067464258930640);
void f_fun(double *state, double dt, double *out_5839489982569705744);
void F_fun(double *state, double dt, double *out_5476922234156819524);
void h_25(double *state, double *unused, double *out_1729658063121940243);
void H_25(double *state, double *unused, double *out_5551778729974662920);
void h_24(double *state, double *unused, double *out_1835918340960958870);
void H_24(double *state, double *unused, double *out_5256232666140554348);
void h_30(double *state, double *unused, double *out_5477911081103537199);
void H_30(double *state, double *unused, double *out_280600472658078466);
void h_26(double *state, double *unused, double *out_5382044204511300316);
void H_26(double *state, double *unused, double *out_1277297089077403520);
void h_27(double *state, double *unused, double *out_4755968270779530176);
void H_27(double *state, double *unused, double *out_5790801706628167388);
void h_29(double *state, double *unused, double *out_2714519565818744062);
void H_29(double *state, double *unused, double *out_5984792589090309848);
void h_28(double *state, double *unused, double *out_3466311643217169430);
void H_28(double *state, double *unused, double *out_8659147327246331634);
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
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
