/* ========================================================================== */
/* == don't modify this part ==*/

/* gsl header files */
//#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>

/* other C header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

/* constants */
#include "constants.h"

/* mnemonics */
#define TRUE (100)
#define FALSE (101)

/* ========================================================================== */
/* == options/parameters related to electron distribution == */

/* mnemonics of distribution function: */
/* THERMAL:
 *      Isotropic relativistic thermal distribution.
 * NONTHERMAL: 
 *      Isotropic nonthermal distribution.
 * NONTHERMAL_ANISO:
 *      Anisotropic nonthermal distribution.
 *      Fleishman & Melnikov 2003 ApJ, 587, 823 Eqs.(11), (13) and (14)
 * GAUSSIAN:
 *      Gaussian function
 *      (1/(DELTA sqrt[Pi])) exp[ - (gamma - GAM0)^2 / DELTA^2 ]
 *      With a small DELTA, can be used to approximate the distribution 
 *      function (Dirac Delta function) of monoenergetic electrons 
 * KAPPA:
 *      Kappa distribution
 *      A deviation of Maxwellian (Thermal) distribution
 *      Kappa distribution converges to Maxwellian distribution at large kappa
 *      (Normalizaiton / (4 * Pi)) * (1 + (sqrt[1 + Pi^2] - 1) / (kappa * theta^2))^(- (kappa + 1))
 *      Normalizaiton term can be determinated both numerically and analytically */
 
#define THERMAL (11)
#define NONTHERMAL (12)
#define NONTHERMAL_ANISO (13)
#define GAUSSIAN (14)
#define KAPPA_DIST (15)
#define KAPPA (150)

/* choose: */
#define ELECTRON_DIST (NONTHERMAL)

/* parameters specific to particular distribution: */
#if (ELECTRON_DIST == THERMAL)
        #define DISTRIBUTION_PARAMETER double ne_para[] = {tdl};
#elif (ELECTRON_DIST == KAPPA_DIST)
        #define DISTRIBUTION_PARAMETER double ne_para[] = {tdl};
#elif (ELECTRON_DIST == NONTHERMAL)
        /* COEFF_INT is Integrate[gamma^(-INDEX), {gamma,GAM1,GAM2}] */
   #if 0
        #define INDEX (3.)
        #define GAM1 (1.)
        #define GAM2 (1000.)
        #define COEFF_INT (0.49999950000000000000)
   #endif
  // #if 0
        #define INDEX (2.2)
        #define GAM1 (1.)
        #define GAM2 (100.)
        #define COEFF_INT (0.8300157735787208562)
   //#endif
   #if 0
        #define INDEX (2.2)
        #define GAM1 (1.)
        #define GAM2 (1000.)
        #define COEFF_INT (0.8331240094640408683)
   #endif
   #if 0
        #define INDEX (2.2)
        #define GAM1 (10.)
        #define GAM2 (100.)
        #define COEFF_INT (0.04926221895207029370)
   #endif
   #if 0
        #define INDEX (2.2)
        #define GAM1 (100.)
        #define GAM2 (10000.)
        #define COEFF_INT (0.003304352311341967811)
   #endif
        #define DISTRIBUTION_PARAMETER double ne_para[] = {INDEX, GAM1, GAM2};
#elif (ELECTRON_DIST == NONTHERMAL_ANISO)
        /* COEFF_INT is Integrate[gamma^(-INDEX), {gamma,GAM1,GAM2}] */
        #define INDEX (5.)
        #define GAM1 (1.01980390271855696600564482180455639791275)
        #define GAM2 (20.02498439450078572769721214832260542149)
        #define COEFF_G_INT (11.59128797788713411539985506775752333303)

        #define XIcFLAG (6)

        /* COEFF_MU_INT is Integrate[f(mu), {mu, -1 , 1}] */
        #if (XIcFLAG == 0)
                #define XIc (0.)
                #define COEFF_MU_INT (2.)
        #elif (XIcFLAG == 2)
                #define XIc (C_pi_2)
                #define COEFF_MU_INT (0.9142857142857142857142857142857142857143)
        #elif (XIcFLAG == 3)
                #define XIc (C_pi/3.)
                #define COEFF_MU_INT (1.484095982142857142857142857142857142857)
        #elif (XIcFLAG == 6)
                #define XIc (C_pi/6.)
                #define COEFF_MU_INT (1.865777494563223052616972008690190423699)
        #else
                #error should not arrive here!
        #endif
        #define DISTRIBUTION_PARAMETER double ne_para[] = {INDEX, GAM1, GAM2, XIc};
#elif (ELECTRON_DIST == GAUSSIAN)
        /* Peak is at gamma = GAM0,
         * FWHM is 1.66511 DELTA */
        #define GAM0 (22.366272042129222)	/* beta = 0.999 -> gamma = 22.366272042129222 */
        //#define GAM0 (1.00727870503172537306)   /* beta = 0.12 -> gamma = 1.00727870503172537306 */
        //#define GAM0 (1.08109687516102642344)   /* beta = 0.38 -> gamma = 1.08109687516102642344 */
        //#define GAM0 (1.9596545041740512626)    /* beta = 0.86 -> gamma = 1.9596545041740512626 */
        //#define GAM0 (1.5)
        #define DELTA (0.1)
        #define DISTRIBUTION_PARAMETER double ne_para[] = {GAM0, DELTA};
#else
        #error Should not arrive here. Invalid value of ELECTRON_DIST.
#endif


/* ========================================================================== */
/* == options/parameters related to angle-averaging emissivity */

/* angle-averaged jnu = (1/2) Integrate[jnu(th) , {cos(th), -1, 1}]
 * Assume that jnu(th) = jnu(-th), therefore
 * angle-averaged jnu = Integrate[jnu(th) , {cos(th), 0, 1}] */

/* number of zone in cos(th) in the integration */
#define N_THETA_AVERAGE (90)

/* ========================================================================== */

/* Mnemonics about different modes
 *  for jnu(), alphanu() and jnu_average() */
#define TOTAL    (0)
#define MODE_O   (1)
#define MODE_X   (-1)
#define STOKES_I (2)
#define STOKES_Q (3)
#define STOKES_V (4)

/* ========================================================================== */
/* == options/parameters related to extrapolation at theta near to Pi/2 */

/* If theta is in (Pi/2 - DELTA_THETA , Pi/2 + DELTA_THETA), then extrapolate.
 *  note that DELTA_THETA is assumed to be positive.  It should also be small
 *  (say, less than 0.0017453 rad or 0.1 deg). */
/* Don't set DELTA_THETA to be too small, or the integrand could be very 
 *  narrow and it would take a lot time to integrate accurately */
#define DELTA_THETA 0.0005
/* Extrapolate result near to Pi/2 ? */
#define EXTRAPOLATE TRUE

/* ========================================================================== */
/* == options/parameters related to integrator == */

/* Mnemonics for choice of integrator */
#define GSL_QAG (0)
#define GSL_QNG (1)

#define INTEGRATOR (GSL_QAG)	/* CHOOSE */

/* note: can't use GSL_QAGIU for n integration although it seems to be a 
 * perfect choice. Turns out GSL_QAGIU can't find a good result because the
 * integrad is not nice and smooth. */

/* -------------------------------------------------------------------------- */
/* Number of double precision intervals held by GSL QAG subroutine */
#define N_WKSP (5000)
//#define N_WKSP (1000)

/* -------------------------------------------------------------------------- */
/* Integration rule used by GSL QAG subroutine */
/* Can be one of GSL_INTEG_GAUSS15, GSL_INTEG_GAUSS21, GSL_INTEG_GAUSS31
 *               GSL_INTEG_GAUSS41, GSL_INTEG_GAUSS51, GSL_INTEG_GAUSS61 */
#define GSL_QAG_RULE_g (GSL_INTEG_GAUSS31)
#define GSL_QAG_RULE_n (GSL_INTEG_GAUSS31)

/* -------------------------------------------------------------------------- */
/* coefficient before Delta n
 * which is used to set the interval of n integral (for thermal distribution)
 * by nEnd = nBegin + COEFF_DELTA_N * some_function */
#define COEFF_DELTA_N 10
//#define COEFF_DELTA_N 40

/* -------------------------------------------------------------------------- */
/* Mnemonics for the way to calculate K_2() */
/* K2_WZ_whole_eq: use the whole fitting formula of K_2 in W&Z 2000,
 * K2_WZ_high_T: use high Theta limit of the formula in W&Z 2000
 * K2_GSL: use GSL
 * K2_GSL_interpolate: use GSL to create a table, then interpolate */
#define K2_WZ_whole_eq (0)
#define K2_WZ_high_T (1)
#define K2_GSL (2)
#define K2_GSL_interpolate (3)

#define K2_EQ (K2_WZ_whole_eq)	/* CHOOSE */
//#define K2_EQ (K2_GSL)	/* CHOOSE */

/* -------------------------------------------------------------------------- */
/* absolute and relative error requirement for gamma and n integration */

/* gamma integral is the integrand of n integral, so gamma integration needs to
 * be more accurate */
#define EPS_REL_n (1.e-6)
//#define EPS_REL_n ((1.e-3*nuratio/1.e5+1.e-8))
#define EPS_REL_g (EPS_REL_n*1.e-2)

/* absolute err should be set to 0 */
#define EPS_ABS_g (0.)
#define EPS_ABS_n (0.)

/* -------------------------------------------------------------------------- */
/* if the distribution is not thermal, intergrate n for
 * [ n_- , n_- + nEND ], where n_- = nu/nuc sin(theta) */
#define nEND (100)

/* In general, use adaptive approach to make sure the integration converges. */
#define ADAPTIVE TRUE

/* parameters of the adaptive approach */
/* NMAX_ADAPTIVE: max number of iteration
 * EPS_ADAPTIVE: max allowed fractional contribution from the last interval */
#define NMAX_ADAPTIVE (200)
//#define EPS_ADAPTIVE (EPS_REL_n*1.e-4)
#define EPS_ADAPTIVE (1.e-8)

/* -------------------------------------------------------------------------- */
/* decide whether use log(n) for independent variable in intergration; 
 * if not, use n */
#define INTEGRATE_WITH_LOGN FALSE

/* -------------------------------------------------------------------------- */
/*  Mnemonics for the way to calculate J_n(x) for small n */
/* JN_C_LIB: use jn(n,x) from c library
 * JN_C_LIB_interpolate: use jn(n,x) from c library, with linear interpolation
 * JN_GSL: use gsl_sf_bessel_Jnu(n,x) from GSL library
 * JN_Chishtie: use Scott's implementation of Chishtie et. al. 2005
 */
#define JN_C_LIB (0)
#define JN_C_LIB_interpolate (1)
#define JN_GSL (2)
#define JN_Chishtie (3)

#define FLAG_N_JN JN_GSL	/* CHOOSE */
//#define FLAG_N_JN JN_C_LIB_interpolate	/* CHOOSE */

/* -------------------------------------------------------------------------- */
/* N_JN: switch from special method (as defined by FLAG_N_JN) to calculate 
 *       Bessel function, to the method by Chishtie et al. (2005) */
#define N_JN (30)
/* n_INT_MIN: minimum n for integration
 *            (should be larger than or equal to N_JN if c library is used
 *            to calculate J_n because n is not continuous in that case */
#define n_INT_MIN (30)

//#if N_SUM < N_JN
//  #error "In decs.h, N_SUM should be larger or equal to N_JN"
//#endif

//#define nterm_SUM_MIN (20)

/* -------------------------------------------------------------------------- */
/* Mnemonics for the way to calculate J'_n(x) */
/* JNprime_EQ1: J'_n = n J_n / x - J_{n+1}
 * JNprime_EQ2: J'_n = 0.5 ( J_{n-1} - J_{n+1} )
 */
#define JNprime_EQ1 (0)
#define JNprime_EQ2 (1)

#define FLAG_JNprime_EQ JNprime_EQ1	/* CHOOSE */ 
/* -------------------------------------------------------------------------- */
/* == d\gamma and d\mu for derivatives in the absorption coefficient 
 *    calculation */
#define DGAMMA (1.e-10)
#define DMU (1.e-10)

/* -------------------------------------------------------------------------- */
/* Limit the boundary of gamma integration for THERMAL distribution, 
 * to exclude unimportant region == */
/* FLAG_GAMMA_BOUNDARY = 0 : use the whole range
 *                       1 : exclude unimportant region */
#define FLAG_GAMMA_BOUNDARY (1)    /* CHOOSE */

/* JN_RATIO : a region is unimportant when 
 *            J_n < JN_RATIO * (maximum value of J_n) */
#define JN_RATIO (1.e-6)

/* sqrt1mr2 : Sqrt[1-r^2] (r is the value of z/n at the new boundary,
 *                         assuming beta=1)
 * nr : number of elements in the array 'sqrt1mr2' */
static double GAMMA_BOUNDARY_sqrt1mr2[] = {0.01,0.05,0.25};
static int GAMMA_BOUNDARY_nr = 3;

/* -------------------------------------------------------------------------- */
/* == Mnemonics for the way to handle unfixable error: == */
/* ERR1_ON_SCREEN_EXIT: print error to stderr and exit
 * ERR1_ON_SCREEN_CONT: print error to stderr and continue
 * ERR1_TO_FILE_CONT  : print error to error.txt and continue
 * ERR1_IGNORE        : ignore error and continue
 */
#define ERR1_ON_SCREEN_EXIT (10)
#define ERR1_ON_SCREEN_CONT (11)
#define ERR1_TO_FILE_CONT   (12)
#define ERR1_IGNORE         (13)

#define FLAG_ERR1 ERR1_IGNORE	/* CHOOSE */

/* == Mnemonics for the way to handle potentially fixable error: == */
/* ERR2_ON_SCREEN_EXIT: print error to stderr and exit
 * ERR2_ON_SCREEN_CONT: print error to stderr, try to fix error, and continue
 * ERR2_TO_FILE_CONT  : print error to error-fixed.txt, try to fix error, and continue
 * ERR2_FIX_CONT      : try to fix error and continue
 * ERR2_IGNORE        : ignore error and continue
 */
#define ERR2_ON_SCREEN_EXIT (20)
#define ERR2_ON_SCREEN_CONT (21)
#define ERR2_TO_FILE_CONT   (22)
#define ERR2_FIX_CONT       (23)
#define ERR2_IGNORE         (24)

#define FLAG_ERR2 ERR2_FIX_CONT	/* CHOOSE */

/*----------------------------------------------------------------------------*/
/* for interpolation of Tx as a function of nu */
typedef struct
{
  int length;
  double * nu;
  double * Tx;
  double rho;
  double Te;
  double B;
  double theta;
  double nuc;
  gsl_interp_accel * acc;
  gsl_spline * spline;
} TYPE_Tx;

/*----------------------------------------------------------------------------*/

typedef struct
{
  double nuratio;
  double tdl;
  double cos_theta;
  double sin_theta;
  double n;
  double B;
  double nrhoe;
  int flag_Tx_interp;
  void * Tx_list;
} TYPE_PAR_g;

/*----------------------------------------------------------------------------*/

typedef struct
{
  TYPE_PAR_g *par_g_pt;
  gsl_function *integrand_g_pt;
  double gBegin;
  double gEnd;
  double gAbsEps;
  double gRelEps;
  size_t n_wksp_g;
  gsl_integration_workspace *wksp_g_pt;
  double gMaxAbsErr;
} TYPE_PAR_n;

/*----------------------------------------------------------------------------*/

typedef struct
{
  double nAbsEps;
  double nRelEps;
  size_t n_wksp_n;
  gsl_integration_workspace *wksp_n_pt;
  double nAbsErr;
  gsl_function integrand_n;
} TYPE_PAR_jnu_I;

/*============================================================================*/
extern double exp_factor(double f_factor, double f_exp);
extern double my_Bessel_dJ(double n, double x, double *jn);
extern double my_Bessel_J(double n, double x);
/*--------------------*/
extern double K2_tdlinv_inv(double tdlinv);
extern double K2_tdlinv_inv_new(double tdlinv);
/*--------------------*/
extern void return_T_L(double *T, double *L, int flag, double nuratio,
                       double nrhoe, double B, double costheta, double sintheta);
/*--------------------*/
extern double jnu(int flagMODE, 
                        double nu, double rho, double Te, double B, double theta,
                        size_t n_wksp, gsl_integration_workspace *wksp_g,
                        gsl_integration_workspace *wksp_n);
/*--------------------*/
double jnu_average(int flagMODE,
                   double nu, double rho, double Te, double B, 
                   size_t n_wksp, gsl_integration_workspace *wksp_g,
                   gsl_integration_workspace *wksp_n);
/*----------------------------------------------------------------------------*/
extern void jnuOX_to_jnu_and_alphanuStokes_Vath
            (double *jnuI, double *jnuQ, double *jnuV,
             double *kappa, double *q, double *v, double *f, double *h,
             double jnuO, double jnuX, double nu, 
             double nrhoe, double Te, double B, double cos_theta, 
             double sin_theta);
extern void jnuOX_to_jnu_and_alphanuStokes_Broderick
            (double *jnuI, double *jnuQ, double *jnuV,
             double *kappa, double *q, double *v,
             double jnuO, double jnuX, double nu, 
             double nrhoe, double Te, double B, double cos_theta, 
             double sin_theta);
extern void jnuOX_to_jnu_and_alphanuStokes_Broderick_new
            (double *jnuI, double *jnuQ, double *jnuV,
             double *kappa, double *q, double *v,
             double jnuO, double jnuX, double nu, 
             double nrhoe, double Te, double B, double cos_theta, 
             double sin_theta);
extern void alphanuOX_to_alphanuStokes_Broderick_new
            (double *kappa, double *q, double *v,
             double alphanuO, double alphanuX, double nuratio, 
             double nrhoe, double B, double cos_theta, double sin_theta);
extern void alphanuStokes_to_alphanuOX_accurate
            (double *alphanuO, double *alphanuX, double *Tx,
             double alphanuI, double alphanuQ, double alphanuV);
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double alphanu(int flagMODE,
                   double nu, double rho, double Te, double B, double theta,
                   size_t n_wksp, gsl_integration_workspace *wksp_g,
                   gsl_integration_workspace *wksp_n,
                   int flag_Tx_interp, void *Tx_interp);
/*============================================================================*/
/* interpolation of Tx: */
void pre_Tx_interp(void *Tx_list, int length, double nu1, double nu2, 
                          double rho, double Te, double B, double theta);
double Tx_interp(void *Tx_list, double nu);
void pro_Tx_interp(void *Tx_list);
/*============================================================================*/

extern double jnu_WZ6(double nu, double nuratio, double nrhoe, double Tdl, double sin_theta);
extern double jnu_T1(double nu, double rho, double Te, double B, double theta);
extern double jnu_T2(double nu, double rho, double Te, double B, double theta);
extern double jnu_T3(double nu, double rho, double Te, double B, double theta);
extern double jnu_T4(double nu, double rho, double Te, double B, double theta);
extern double jnu_T5(double nu, double rho, double Te, double B, double theta);
extern double jnu_T6(double nu, double rho, double Te, double B, double theta);

/* vim:set ts=8 expandtab: */
