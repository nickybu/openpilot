/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7455540085193669325);
void inv_err_fun(double *nom_x, double *true_x, double *out_918926006996616643);
void H_mod_fun(double *state, double *out_4843743724872727942);
void f_fun(double *state, double dt, double *out_4710965010460499451);
void F_fun(double *state, double dt, double *out_2935311077495091398);
void h_3(double *state, double *unused, double *out_5619078680604814377);
void H_3(double *state, double *unused, double *out_1189855813453559207);
void h_4(double *state, double *unused, double *out_7052672877425575945);
void H_4(double *state, double *unused, double *out_8975370522166609105);
void h_9(double *state, double *unused, double *out_2759646287952145907);
void H_9(double *state, double *unused, double *out_4889524503794664352);
void h_10(double *state, double *unused, double *out_4001570132546139827);
void H_10(double *state, double *unused, double *out_7828488975012777857);
void h_12(double *state, double *unused, double *out_515585852484766857);
void H_12(double *state, double *unused, double *out_6108065045044855600);
void h_31(double *state, double *unused, double *out_278764678192433924);
void H_31(double *state, double *unused, double *out_2689252026114668008);
void h_32(double *state, double *unused, double *out_5164909157540339636);
void H_32(double *state, double *unused, double *out_3877789004048678586);
void h_13(double *state, double *unused, double *out_4359963420172572358);
void H_13(double *state, double *unused, double *out_1286944411026122826);
void h_14(double *state, double *unused, double *out_2759646287952145907);
void H_14(double *state, double *unused, double *out_4889524503794664352);
void h_19(double *state, double *unused, double *out_7442283347674055381);
void H_19(double *state, double *unused, double *out_6921733424903324572);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);