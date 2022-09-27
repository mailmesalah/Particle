/*****************************
  PARTICLE MODEL 1-D PYROLYSIS 
  TU BERLIN - ANDRES ANCA COUCE
 *****************************/

/************************
  DECLARATION OF VARIABLES
 ************************/
#define NGRID 20
#define NGSP 5
#define NSSP 2
#define NRXN 4
#define NEQ 8      
/* NGRID = number of finite volumes in the particle
   NGSP = number of gas species
   NSSP = number of solid species
   NRXN = number of reactions
   NEQ = number of ODE in each finite volume */

int OPT_1, OPT_2;
/* OPT_2 = Difussion, 1:Fick, 2:Dusty gas model, 3:Wilke, else:no
   OPT_1 = Geometrical configuration, 1:spherical, else:Gronli */

double L, RAD, H0;
double SHAPE_IN[NGRID+1], SHAPE_OUT[NGRID+1];
/* L = Lenght of the longitudinal particle
   RAD = Radius of the spherical particle
   H0 = Lenght of the finite volume
   SHAPE_IN = Spherical shape factor of the inside boundary  
   SHAPE_OUT = Spherical shape factor of the outside boundary
   Shape factor is the relation between the area of the surface 
   and the volume of the finite element */

/* SHR_FACTOR_MIN = Minimum shrinkage factor (-)
   DENSITY = Density of the gas mixture (kg/m3)
   PERMEABILITY = Intinsic permeability of the solid (m2)
   D_POR = Diameter of the por (m)
   THERM_CONDUC = Thermal conductivity (W/m*K)
   EMISS = Solid emissivity (-)
   CP = Specific heat capacity (J/kg*K)
   WM = Molecular weight (kg/mol)
   VISCOSITY = Dynamic viscosity of the gas mixture (kg/m*s)
   DIFUSS_FACTOR = Factor of porous structure in difussion (-)
   STF_BOL = Stefan-Boltzman constant (W/m2*K4)
   CTE_GAS = Gas constant (J/K*mol) */
      
double SHR_FACTOR_MIN;
double DENSITY_SOLID[NSSP+1];
double PERMEABILITY_WOOD, PERMEABILITY_CHAR, D_POR_WOOD, D_POR_CHAR;
double THERM_CONDUC_WOOD, THERM_CONDUC_CHAR, EMISS, STF_BOL;
double CP_SSP_0[NSSP+1], CP_SSP_1[NSSP+1];
double WM_GSP[NGSP+1], CP_GSP_0[NGSP+1], CP_GSP_1[NGSP+1], CP_GSP_2[NGSP+1];
double D_SIGMA[NGSP+1], D_EPSILON_K[NGSP+1];
double D_AB[NGSP+1][NGSP+1];
double VISCOSITY_GAS, THERM_CONDUC_GAS, DIFUSS_FACTOR, CTE_GAS;
double Y_GAS_OUT[NGSP+1], X_GAS_OUT[NGSP+1];
double T_OUT, PRES_OUT, WM_GAS_OUT, DENSITY_GAS_OUT;
double F_FLUX, ALPHA, BETA, S_EMISS;
double POROSITY_0, T_0, P_0, WM_GAS_0, DENSITY_GAS_0;
double Y_GAS_0[NGSP+1], Y_SOLID_0[NSSP+1];
double SCALE_SOLID, SCALE_T;

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double getTime();
void inputs(void);

int main(void) {
  long int I, SPC, R, S, DUM;
  double DUMMY;
  double CONV[NGRID+1], SHR_FACTOR[NGRID+1], RW[NGRID+1];
  double H_P[NGRID+1], RHO_POROS_SOLID[NGRID+1], POROS[NGRID+1];
  double RHO_POROS_GAS[NGRID+1], DENSITY_GAS[NGRID+1];
  double WM_GAS[NGRID+1], T[NGRID+1], PRES[NGRID+1];
  double CP_GSP[NGSP+1], CP_SSP[NSSP+1];
  double CP_GAS[NGRID+1],RHO_CP[NGRID+1];
  double PERMEABILITY[NGRID+1], CONDUC_TOT[NGRID+1];
  double Y_GAS[NGRID+1][NGSP+1], X_GAS[NGRID+1][NGSP+1];
  double D_POR, CONDUC_RAD;
  double PERMEABILITY_S, CONDUC_S, T_S_EXTRAP, T_S;
  double EPSILON_K, T_PRIMA, OMEGA;
  double DIFF_AB[NGRID+1][NGSP+1][NGSP+1];
  double DIFF[NGRID+1][NGSP+1];
  double RATE[NGRID+1][NRXN+1], R_GAS_TOT[NGRID+1];
  double R_GAS[NGRID+1][NGSP+1], R_SOLID[NGRID+1][NSSP+1];
  double RXN_ENTH[NRXN+1];
  double HEAT_REACTION;

  double TNOW, T_STOP, DELTA_T;
  double Y_VECTOR[NGRID+1][NEQ+1], Y_VECTOR_BTS[NGRID+1][NEQ+1];
  long int N_STEP, MAX_N_STEP, N_ITER, MAX_N_ITER;
  double SUM_CORREC, MIN_SUM_CORREC;

  double DENSITY_GAS_BTS[NGRID+1], WM_GAS_BTS[NGRID+1];
  double PRES_BTS[NGRID+1], POROS_BTS[NGRID+1], T_BTS[NGRID+1];

  double A_IN, A_OUT, CONDUC_IN, CONDUC_OUT;
  double AP[NGRID+1], AW[NGRID+1], AE[NGRID+1],Q[NGRID+1];
  double PRES_PRED[NGRID+1], VEL[NGRID+1], VEL_C[NGRID+1];

  double TDMA[NGRID+1];
  double aa[NGRID+1], bb[NGRID+1], cc[NGRID+1], dd[NGRID+1];
  double start, finish;

  inputs();

  // Capture the starting time  
  start = getTime();

  T_STOP = 100.0;
  DELTA_T = 0.0001;
  MAX_N_STEP = 10000;

  TNOW = 0.0;

  for (I = 1; I <= NGRID; I++) {
    for (SPC = 1; SPC <= NGSP; SPC++) {
      Y_VECTOR[I][SPC] = POROSITY_0 * DENSITY_GAS_0 * Y_GAS_0[SPC];
    }
    for (SPC = 1; SPC <= NSSP; SPC++) {
      Y_VECTOR[I][NGSP+SPC] = (1.0 - POROSITY_0) * DENSITY_SOLID[SPC] * Y_SOLID_0[SPC];
    }
    Y_VECTOR[I][NGSP+NSSP+1] = T_0;
  }

  #include "eq_state.h"

/******************************************
 --------------- SIMULATION ---------------
 ******************************************/

  for (N_STEP = 1; N_STEP <= MAX_N_STEP; N_STEP++) {
    TNOW = TNOW + DELTA_T;

/*****************
  PROPERTIES UPDATE
 *****************/

    for (I = 1; I <= NGRID; I++) {
          
/* Conversion and shrinkage factors */
      CONV[I] = Y_VECTOR[I][NGSP+1] / ((1.0 - POROSITY_0) * DENSITY_SOLID[1]);
      SHR_FACTOR[I] = SHR_FACTOR_MIN + CONV[I] * (1.0 - SHR_FACTOR_MIN);
         
/*  Lenght of the control volume and shape factors */
      if (OPT_1 == 1) {
        double ARG;
        if (I == 1) {
          RW[I] = 0.0;
        } else {
          RW[I] = RW[I-1] + H_P[I-1];
        }
        ARG = H0*H0*H0 * (3.0 * I*I - 3.0*I + 1.0);
        H_P[I] = pow((RW[I]*RW[I]*RW[I] + SHR_FACTOR[I]*ARG),(1.0/3.0)) - RW[I];
        SHAPE_IN[I] = 3.0 * RW[I]*RW[I] / ARG;
        SHAPE_OUT[I] = 3.0 * pow((RW[I] + H_P[I]),2.0) / ARG;
      } else {
        H_P[I] = H0 * SHR_FACTOR[I];
        SHAPE_IN[I] = 1.0 / H0;
        SHAPE_OUT[I] = 1.0 / H0;
      }
    }
          
    #include "eq_state.h"

/* Get values beginning of the time step */

    for (I=1; I<=NGRID; I++) {
      for (SPC=1; SPC<=NEQ; SPC++) {
        Y_VECTOR_BTS[I][SPC] = Y_VECTOR[I][SPC];
      }
      T_BTS[I] = T[I];
      PRES_BTS[I] = PRES[I];
      DENSITY_GAS_BTS[I] = DENSITY_GAS[I];
      WM_GAS_BTS[I] = WM_GAS[I];
      POROS_BTS[I] = POROS[I];
    }

    for (I=1; I<=NGRID; I++) {

//  Heat capacity
      for (SPC=1; SPC<=NGSP; SPC++) {
        CP_GSP[SPC] = CP_GSP_0[SPC] + CP_GSP_1[SPC] * T[I] + CP_GSP_2[SPC] * T[I]*T[I];
      } 
      for (SPC=1; SPC<=NSSP; SPC++) {
        CP_SSP[SPC] = CP_SSP_0[SPC] + CP_SSP_1[SPC] * T[I];
      }
      RHO_CP[I] = 0.0;
      for (SPC=1; SPC<=NGSP; SPC++) {
        RHO_CP[I] = RHO_CP[I] + CP_GSP[SPC] * Y_VECTOR[I][SPC];
      }
      CP_GAS[I] = RHO_CP[I] / (POROS[I] * DENSITY_GAS[I]);
      for (SPC=1; SPC<=NSSP; SPC++) {
        RHO_CP[I] = RHO_CP[I] + CP_SSP[SPC] * Y_VECTOR[I][NGSP+SPC];
      }

// Permeability
      PERMEABILITY[I] = CONV[I] * PERMEABILITY_WOOD + (1.0 - CONV[I]) * PERMEABILITY_CHAR;
    
// Thermal conductivity
      D_POR = CONV[I] * D_POR_WOOD + (1.0 - CONV[I]) * D_POR_CHAR;
      CONDUC_RAD = 4.0 * POROS[I] * STF_BOL * EMISS * D_POR \
            * (T[I]*T[I]*T[I]) / ( 1.0 - POROS[I]);
     
      CONDUC_TOT[I] = (POROS[I] * THERM_CONDUC_GAS \
            + CONV[I] * THERM_CONDUC_WOOD \
            + (1 - CONV[I]) * THERM_CONDUC_CHAR + CONDUC_RAD);

    }
      
/*------------------
  Surface properties
  ------------------*/

    I = NGRID;
      
//  Extrapolated surface properties
    PERMEABILITY_S = 1.5 * PERMEABILITY[NGRID] - 0.5 * PERMEABILITY[NGRID-1];
    CONDUC_S = 1.5 * CONDUC_TOT[NGRID] - 0.5 * CONDUC_TOT[NGRID-1];
    T_S_EXTRAP = 1.5 * T[NGRID] - 0.5 * T[NGRID-1];
      
// Energy balace in the surface
    T_S = ( F_FLUX + ALPHA * T_OUT + S_EMISS * STF_BOL * pow(T_OUT,4.0) \
        + 2.0 * CONDUC_TOT[NGRID] * T[NGRID] / H_P[NGRID] ) \
        / ( ALPHA + S_EMISS * STF_BOL * pow(T[NGRID],3.0) \
        + 2.0 * CONDUC_TOT[NGRID] / H_P[NGRID] );

/***************
  DIFUSSIVE TERMS
 ***************/
    
// Matrix binary diffusion coefficients DIFF_AB, calculated surface between CV
// border i between CV i and i+1

    for (I=1; I<=(NGRID-1); I++) {
      for (R=1; R<=NGSP; R++) {
        for (S=1; S<=(R-1); S++) {
          EPSILON_K = sqrt(D_EPSILON_K[R] * D_EPSILON_K[S]);
          T_PRIMA = (T[I] + T[I+1]) / (2.0 * EPSILON_K);
// Omega: diffusion collision integral, dimensionless
          OMEGA = 1.06036 / pow(T_PRIMA,0.15610) \
                  + 0.193 / exp(0.47635*T_PRIMA) \
                  + 1.03587 / exp(1.52996*T_PRIMA) \
                  + 1.76474 / exp(3.89411*T_PRIMA);
/* D_AB: diffusion coefficient, m2/s
   Considering pressure: 1 bar; include tortuosity (DIFFUS_FACTOR)
   D_AB(i,j) = 1.0e-4 * 0.00266 * T^(3.0/2.0) / ...
   (1.0 * WM_AB^0.5 * SIGMA_AB^2.0 * OMEGA); */
          DIFF_AB[I][R][S] = \
          D_AB[R][S] * pow((T[I] + T[I+1])/2.0,(3.0/2.0)) / (OMEGA * DIFUSS_FACTOR);
// D_AB is triangular inferior matrix, not DIFF_AB
          DIFF_AB[I][S][R] = DIFF_AB[I][R][S];
        }
        // Knudsen diffusion
        DIFF_AB[I][R][R] = 0.0; //Check it
      }
    }

    I = NGRID;
    for (R=1; R<=NGSP; R++) {
      for (S=1; S<=(R-1); S++) {
        EPSILON_K = sqrt(D_EPSILON_K[R] * D_EPSILON_K[S]);
        T_PRIMA = T_S / (2.0 * EPSILON_K);
    
        OMEGA = 1.06036 / pow(T_PRIMA,0.15610) \
                + 0.193 / exp(0.47635*T_PRIMA) \
                + 1.03587 / exp(1.52996*T_PRIMA) \
                + 1.76474 / exp(3.89411*T_PRIMA);
  
        DIFF_AB[I][R][S] = D_AB[R][S] * pow(T_S,(3.0/2.0)) / (OMEGA * DIFUSS_FACTOR);
        DIFF_AB[I][S][R] = DIFF_AB[I][R][S];
      }     
      // Knudsen diffusion
      DIFF_AB[I][R][R] = 0.0; //Check it
    }

/*  Mass diffusion flux: J_DIFF(kg/m2*s)
    Calculated in the border between CV, border i between CV i and i+1
    No flux in simmetry axis */

    if (OPT_2 == 1) {
// *** Fickian diffusion ***
      for (I=1; I<=NGRID; I++) {
        for (SPC=1; SPC<=(NGSP-1); SPC++) {
          DIFF[I][SPC] = DIFF_AB[I][SPC][NGSP];
        }
      }
    }

/********
  REACTION
 ********/

    for (I=1; I<=NGRID; I++) {
//  RATE = Reaction rate 1 (kg/m3*s)
      RATE[I][1] = Y_VECTOR[I][NGSP+1] * 1.3e8 * exp(-140.3e3 / (CTE_GAS * T[I]));
      RATE[I][2] = Y_VECTOR[I][NGSP+1] * 2.0e8 * exp(-133.1e3 / (CTE_GAS * T[I]));
      RATE[I][3] = Y_VECTOR[I][NGSP+1] * 1.1e7 * exp(-121.3e3 / (CTE_GAS * T[I]));
      RATE[I][4] = Y_VECTOR[I][2] * 4.3e6 * exp(-108.0e3 / (CTE_GAS * T[I]));
      R_GAS_TOT[I] = RATE[I][1] + RATE[I][2];
    }

/*--------------------------------
  ---------- ITERATIONS ----------
  --------------------------------*/

    MAX_N_ITER = 20;
    MIN_SUM_CORREC = 0.01;
    SUM_CORREC = 1.0;

// while ( N_ITER < MAX_N_ITER & SUM_CORREC > MIN_SUM_CORREC )
		
    for (N_ITER=1; N_ITER<=MAX_N_ITER; N_ITER++) {
		
/*****************
  PRESSURE EQUATION
 *****************/

// Predicted pressure

      I = 1;

      A_IN = 0.0;
      if (PRES[I] >= PRES[I+1]) {    
        A_OUT = SHAPE_OUT[I] * DENSITY_GAS_BTS[I] * \
             ((PERMEABILITY[I] + PERMEABILITY[I+1]) / 2.0) / \
             ((H_P[I]+H_P[I+1]) * VISCOSITY_GAS / 2.0);
      } else {
        A_OUT = SHAPE_OUT[I] * DENSITY_GAS_BTS[I+1] *
             ((PERMEABILITY[I] + PERMEABILITY[I+1]) / 2.0) / \
             ((H_P[I]+H_P[I+1]) * VISCOSITY_GAS / 2.0);
      }
    
      AP[I] = (SHR_FACTOR[I] * POROS[I] * WM_GAS[I]) /  \
              (CTE_GAS * T[I] * DELTA_T) + A_IN + A_OUT;
      AW[I] = -A_IN;
      AE[I] = -A_OUT;
      Q[I] = (SHR_FACTOR[I] * POROS_BTS[I]*WM_GAS_BTS[I]*PRES_BTS[I]) / \
             (CTE_GAS * T_BTS[I] * DELTA_T) + R_GAS_TOT[I];

      for (I=2; I<=(NGRID-1); I++) {
         
        if (PRES[I-1] >= PRES[I]) {    
          A_IN = SHAPE_IN[I] * DENSITY_GAS_BTS[I-1] * \
                 ((PERMEABILITY[I-1] + PERMEABILITY[I]) / 2.0) / \
                 ((H_P[I-1]+H_P[I]) * VISCOSITY_GAS / 2.0);
        } else {
          A_IN = SHAPE_IN[I] * DENSITY_GAS_BTS[I] * \
                 ((PERMEABILITY[I-1] + PERMEABILITY[I]) / 2.0) / \
                 ((H_P[I-1]+H_P[I]) * VISCOSITY_GAS / 2.0);
        }
         
        if (PRES[I] >= PRES[I+1]) {    
          A_OUT = SHAPE_OUT[I] * DENSITY_GAS_BTS[I] * \
                 ((PERMEABILITY[I] + PERMEABILITY[I+1]) / 2.0) / \
                 ((H_P[I]+H_P[I+1]) * VISCOSITY_GAS / 2.0);
        } else {
          A_OUT = SHAPE_OUT[I] * DENSITY_GAS_BTS[I+1] * \
                 ((PERMEABILITY[I] + PERMEABILITY[I+1]) / 2.0) / \
                 ((H_P[I]+H_P[I+1]) * VISCOSITY_GAS / 2.0);
        }
         
        AP[I] = (SHR_FACTOR[I] * POROS[I] * WM_GAS[I]) / \
                (CTE_GAS * T[I] * DELTA_T) + A_IN + A_OUT;
        AW[I] = -A_IN;
        AE[I] = -A_OUT;
        Q[I] = (SHR_FACTOR[I]*POROS_BTS[I]*WM_GAS_BTS[I]*PRES_BTS[I]) / \
               (CTE_GAS * T_BTS[I] * DELTA_T) + R_GAS_TOT[I];
      }
      
      I = NGRID;
      
      if ( PRES[I-1] >= PRES[I]) {    
        A_IN = SHAPE_IN[I] * DENSITY_GAS_BTS[I-1] * \
              ((PERMEABILITY[I-1] + PERMEABILITY[I]) / 2.0) / \
              ((H_P[I-1]+H_P[I]) * VISCOSITY_GAS / 2.0);
      } else {
        A_IN = SHAPE_IN[I] * DENSITY_GAS_BTS[I] * \
              ((PERMEABILITY[I-1] + PERMEABILITY[I]) / 2.0) / \
              ((H_P[I-1]+H_P[I]) * VISCOSITY_GAS / 2.0);
      }
      
      if (PRES[I] >= PRES_OUT) {
        A_OUT = SHAPE_OUT[I] * DENSITY_GAS_BTS[I] * PERMEABILITY_S * 2.0 / \
                (VISCOSITY_GAS * H_P[NGRID]);
      } else {
        A_OUT = SHAPE_OUT[I] * DENSITY_GAS_OUT * PERMEABILITY_S * 2.0 / \
                (VISCOSITY_GAS * H_P[NGRID]);
      }
    
      AP[I] = (SHR_FACTOR[I] * POROS[I] * WM_GAS[I]) / \
              (CTE_GAS * T[I] * DELTA_T) + A_IN + A_OUT;
      AW[I] = -A_IN;
      AE[I] = 0.0;
      Q[I] = (SHR_FACTOR[I] * POROS_BTS[I] * WM_GAS_BTS[I] * PRES_BTS[I]) / \
             (CTE_GAS*T_BTS[I]*DELTA_T) + R_GAS_TOT[I] + A_OUT * PRES_OUT;

// TDMA solver
      for (I=1; I<=NGRID; I++) {
         aa[I] = AW[I];
         bb[I] = AP[I];
         cc[I] = AE[I];
         dd[I] = Q[I];
      }
      I = 1;
      cc[I] = cc[I] / bb[I];
      dd[I] = dd[I] / bb[I];
      for (I=2; I<=NGRID; I++) {
         cc[I] = cc[I] / (bb[I] - cc[I-1]*aa[I]);
         dd[I] = (dd[I] - dd[I-1]*aa[I]) / (bb[I] - cc[I-1]*aa[I]);
      }
      I = NGRID;
      TDMA[I] = dd[I];
      for (I=(NGRID-1); I>=1; I--) {
        TDMA[I] = dd[I] - cc[I] * TDMA[I+1];
      }
      
      for (I=1; I<=NGRID; I++) {
        PRES_PRED[I] = TDMA[I];
      }
	  
      I = 1;
      VEL[I] = 0.0;
      
      for (I=2; I<=NGRID; I++) {
        VEL[I] = -((PERMEABILITY[I] + PERMEABILITY[I-1]) / 2.0) * \
                   (PRES_PRED[I] - PRES_PRED[I-1]) / \
                   ((H_P[I-1]+H_P[I]) * VISCOSITY_GAS / 2.0);
      }

      I = NGRID+1;
      VEL[I] = -PERMEABILITY_S * (PRES_OUT - PRES_PRED[I-1]) * 2.0 / (H_P[NGRID] * VISCOSITY_GAS);
           
      for (I=1; I<=NGRID; I++) {   
        VEL_C[I] = (VEL[I] + VEL[I+1]) / 2.0;
      }

/**********************************
  INTEGRATION MASS SPECIES EQUATIONS
  **********************************/

/*----------
  CONVECTION
  ----------*/
// Explicit convection

      I = 1;
      for (SPC=1; SPC<=NGSP; SPC++) {
        Y_VECTOR[I][SPC] = Y_VECTOR_BTS[I][SPC] + DELTA_T * \
             (-SHAPE_OUT[I] * max(VEL[I+1],0.0) * Y_VECTOR_BTS[I][SPC] \
             / (POROS_BTS[I] * SHR_FACTOR[I]) \
             + SHAPE_OUT[I] * max(-VEL[I+1],0.0) * Y_VECTOR_BTS[I+1][SPC] \
             / (POROS_BTS[I+1] * SHR_FACTOR[I+1]));
      }
      
      for (I=2; I<=(NGRID-1); I++) {
         for (SPC=1; SPC<=NGSP; SPC++) {
           Y_VECTOR[I][SPC] = Y_VECTOR_BTS[I][SPC] + DELTA_T * 
                 (SHAPE_IN[I] * max(VEL[I],0.0) * Y_VECTOR_BTS[I-1][SPC] \
                 / (POROS_BTS[I-1] * SHR_FACTOR[I-1]) \
                 - SHAPE_IN[I] * max(-VEL[I],0.0) * Y_VECTOR_BTS[I][SPC] \
                 / (POROS_BTS[I] * SHR_FACTOR[I]) \
                 - SHAPE_OUT[I] * max(VEL[I+1],0.0) * Y_VECTOR_BTS[I][SPC] \
                 / (POROS_BTS[I] * SHR_FACTOR[I]) \
                 + SHAPE_OUT[I] * max(-VEL[I+1],0.0) * Y_VECTOR_BTS[I+1][SPC] \
                 / (POROS_BTS[I+1] * SHR_FACTOR[I+1]));
        }
      }

      I = NGRID;
      for (SPC=1; SPC<=NGSP; SPC++) {
         Y_VECTOR[I][SPC] = Y_VECTOR_BTS[I][SPC] + DELTA_T * \
               (SHAPE_IN[I] * max(VEL[I],0.0) * Y_VECTOR_BTS[I-1][SPC] \
             / (POROS_BTS[I-1] * SHR_FACTOR[I-1]) \
             - SHAPE_IN[I] * max(-VEL[I],0.0) * Y_VECTOR_BTS[I][SPC] \
             / (POROS_BTS[I] * SHR_FACTOR[I]) \
             - SHAPE_OUT[I] * max(VEL[I+1],0.0) * Y_VECTOR_BTS[I][SPC] \
             / (POROS_BTS[I] * SHR_FACTOR[I]) \
             + SHAPE_OUT[I] * max(-VEL[I+1],0.0) \
             * DENSITY_GAS_OUT * Y_GAS_OUT[SPC]);
      }
  
/*--------
  REACTION
  --------*/
// Explicit pyrolysis
     
      for (I=1; I<=NGRID; I++) {
/* R_GAS = Rate of formation of other gases (kg/m3*s)
   Gas pyrolysis; 1:Other gas, 2:Tar1, 3:Tar2, 4:H2O, 5:Dry air  */
        R_GAS[I][1] = RATE[I][1] + RATE[I][4];
        R_GAS[I][2] = RATE[I][2] - RATE[I][4];
        R_GAS[I][3] = 0.0;
        R_GAS[I][4] = 0.0;
        R_GAS[I][5] = 0.0;
        R_SOLID[I][1] = -RATE[I][1] - RATE[I][2] - RATE[I][3];
        R_SOLID[I][2] = RATE[I][3];
      }
      
/* DELTAH = Heat of reaction (J/kg) */
      
      RXN_ENTH[1] = 150.0e3;
      RXN_ENTH[2] = 150.0e3;
      RXN_ENTH[3] = 150.0e3;
      RXN_ENTH[4] = -50.0e3;

      for (I=1; I<=NGRID; I++) {
        for (SPC=1; SPC<=NGSP; SPC++) {
          Y_VECTOR[I][SPC] = Y_VECTOR[I][SPC] + DELTA_T * R_GAS[I][SPC];
        }
        for (SPC=1; SPC<=NSSP; SPC++) {
          Y_VECTOR[I][NGSP+SPC] = Y_VECTOR_BTS[I][NGSP+SPC] + DELTA_T * R_SOLID[I][SPC];
        }
      }

/***************************
  INTEGRATION ENERGY EQUATION
 ***************************/

/*---------
  ADVECTION
  ---------*/
// Explicit advection

      I = 1;
      Y_VECTOR[I][NEQ] = Y_VECTOR_BTS[I][NEQ] + \
        (DELTA_T * DENSITY_GAS_BTS[I] * CP_GAS[I] / RHO_CP[I]) * \
        (-SHAPE_OUT[I] * max(VEL_C[I],0.0) * Y_VECTOR_BTS[I][NEQ] \
        + SHAPE_OUT[I] * max(-VEL_C[I],0.0) * Y_VECTOR_BTS[I+1][NEQ]);

      for (I=2; I<=(NGRID-1); I++) {   
        Y_VECTOR[I][NEQ] = Y_VECTOR_BTS[I][NEQ] + \
             (DELTA_T * DENSITY_GAS_BTS[I] * CP_GAS[I] / RHO_CP[I]) * \
             ( SHAPE_IN[I] * max(VEL_C[I],0.0) * Y_VECTOR_BTS[I-1][NEQ] \
             - SHAPE_IN[I] * max(-VEL_C[I],0.0) * Y_VECTOR_BTS[I][NEQ] \
             - SHAPE_OUT[I] * max(VEL_C[I],0.0) * Y_VECTOR_BTS[I][NEQ] \
             + SHAPE_OUT[I] * max(-VEL_C[I],0.0) * Y_VECTOR_BTS[I+1][NEQ]);  
      }

      I = NGRID;
      Y_VECTOR[I][NEQ] = Y_VECTOR_BTS[I][NEQ] + 
          ( DELTA_T * DENSITY_GAS_BTS[I] * CP_GAS[I] / RHO_CP[I] ) *
          ( SHAPE_IN[I] * max(VEL_C[I],0.0) * Y_VECTOR_BTS[I-1][NEQ] \
          - SHAPE_IN[I] * max(-VEL_C[I],0.0) * Y_VECTOR_BTS[I][NEQ] \
          - SHAPE_OUT[I] * max(VEL_C[I],0.0) * Y_VECTOR_BTS[I][NEQ] \
          + SHAPE_OUT[I] * max(-VEL_C[I],0.0) * T_OUT);

/*----------
  CONDUCTION
  ----------*/
// Conduction semi-implicit
      
      I = 1;
      CONDUC_IN = 0.0;
      CONDUC_OUT = SHAPE_OUT[I] * ((CONDUC_TOT[I] + CONDUC_TOT[I+1]) / 2.0) / \
          ((H_P[I]+H_P[I+1]) / 2.0);
    
      AP[I] = RHO_CP[I] / DELTA_T + CONDUC_IN + CONDUC_OUT;
      AW[I] = -CONDUC_IN;
      AE[I] = -CONDUC_OUT;
      Q[I] = RHO_CP[I] * Y_VECTOR[I][NEQ] / DELTA_T;

      for (I=2; I<=(NGRID-1); I++) {
        CONDUC_IN = SHAPE_IN[I] * ((CONDUC_TOT[I-1] + CONDUC_TOT[I]) / 2.0) / \
          ((H_P[I-1]+H_P[I]) / 2.0);
        CONDUC_OUT = SHAPE_OUT[I] * ((CONDUC_TOT[I] + CONDUC_TOT[I+1]) / 2.0) / \
          ((H_P[I]+H_P[I+1]) / 2.0);
    
         AP[I] = RHO_CP[I] / DELTA_T + CONDUC_IN + CONDUC_OUT;
         AW[I] = -CONDUC_IN;
         AE[I] = -CONDUC_OUT;
         Q[I] = RHO_CP[I] * Y_VECTOR[I][NEQ] / DELTA_T;
      }

      I = NGRID;
      CONDUC_IN = SHAPE_IN[I] * ((CONDUC_TOT[I-1] + CONDUC_TOT[I]) / 2.0) / \
         ((H_P[I-1]+H_P[I]) / 2.0);
      CONDUC_OUT = SHAPE_OUT[I] * CONDUC_S * 2.0 / H_P[NGRID];
    
      AP[I] = RHO_CP[I] / DELTA_T + CONDUC_IN + CONDUC_OUT;
      AW[I] = -CONDUC_IN;
      AE[I] = 0.0;
      Q[I] = RHO_CP[I] * Y_VECTOR[I][NEQ] / DELTA_T + CONDUC_OUT * T_S;

// TDMA solver
      for (I=1; I<=NGRID; I++) {
         aa[I] = AW[I];
         bb[I] = AP[I];
         cc[I] = AE[I];
         dd[I] = Q[I];
      }
      I = 1;
      cc[I] = cc[I] / bb[I];
      dd[I] = dd[I] / bb[I];
      for (I=2; I<=NGRID; I++) {
        cc[I] = cc[I] / (bb[I] - cc[I-1]*aa[I]);
        dd[I] = (dd[I] - dd[I-1]*aa[I]) / (bb[I] - cc[I-1]*aa[I]);
      }
      I = NGRID;
      TDMA[I] = dd[I];
      for (I=(NGRID-1); I>=1; I--) {
         TDMA[I] = dd[I] - cc[I] * TDMA[I+1];
      }
      
      for (I=1; I<=NGRID; I++) {
         Y_VECTOR[I][NEQ] = TDMA[I];
      }

/*--------
  REACTION
  --------*/

      for (I=1; I<=NGRID; I++) {
         HEAT_REACTION = 0.0;
         for (SPC=1; SPC<=NRXN; SPC++) {
            HEAT_REACTION = HEAT_REACTION + RATE[I][SPC] * RXN_ENTH[SPC];
         }
         Y_VECTOR[I][NEQ] = Y_VECTOR[I][NEQ] - ( DELTA_T / RHO_CP[I]) * HEAT_REACTION;
      }

/*-----------------
  CHECK CONVERGENCE
  ----------------*/

      #include "eq_state.h"

      SUM_CORREC = 0.0;
      for (I=1; I<=NGRID; I++) {
         SUM_CORREC = SUM_CORREC + fabs(PRES[I]-PRES_PRED[I]);
      }

      if (SUM_CORREC <= MIN_SUM_CORREC) {
        goto end_iterations;
      }

    }

/*--------------
  End iterations
  --------------*/
end_iterations:
    if (N_ITER >= MAX_N_ITER) {
      printf("Iterations limit reached at time: %lf \t %lf\n", TNOW,SUM_CORREC);
      exit(1);
    }

  }

  finish = getTime();
  printf("Time step (s): %lf\n", DELTA_T);
  printf("Execution time (s): %lf\n", finish-start);

  printf("%d\n", N_ITER);

  for (I=1; I<=NGRID; I++) {
    printf("%.15lf\n", Y_VECTOR[I][NGSP+NSSP+1]);
  }
  for (I=1; I<=NGRID; I++) {
    printf("%.15lf\n", Y_VECTOR[I][1]);
  }

  return 0;
}

void inputs() {
  int I, SPC;
  double WM_AB, SIGMA_AB;

  OPT_1 = 1;
  OPT_2 = 0;
      
  if (OPT_1 == 1) {
/* SPHERICAL CONFIGURATION */
    RAD = 0.010;
    H0 = RAD / NGRID;
  } else {
/* LONGITUDINAL CONFIGURATION */
    L = 0.01;
    H0 = L / NGRID;
  }
      
  SHR_FACTOR_MIN = 0.5;

/* Solid; 1:Wood, 2:Char */
  DENSITY_SOLID[1] = 1400.0;
  DENSITY_SOLID[2] = 1540.0;

  PERMEABILITY_WOOD = 1.0e-14;
  PERMEABILITY_CHAR = 1.0e-12;

  D_POR_WOOD = 50.0e-6;
  D_POR_CHAR = 100.0e-6;

  THERM_CONDUC_WOOD = 0.25;
  THERM_CONDUC_CHAR = 0.1;

  STF_BOL = 5.67e-8;
  EMISS = 0.9;
	  
/* The heat capacity changes with temperature
  CP_SSP(j) = CP_SSP_0(j) + CP_SSP_1(j) * T[I] */
  CP_SSP_0[1] = 1500.0;
  CP_SSP_1[1] = 1.0;
	   
  CP_SSP_0[2] = 808.9; 
  CP_SSP_1[2] = 0.93;

/* Gas pyrolysis; 1:Other gas, 2:Tar1, 3:Tar2, 4:H2O, 5:Dry air */
  WM_GSP[1] = 36e-3;
  WM_GSP[2] = 110e-3;
  WM_GSP[3] = 110e-3;
  WM_GSP[4] = 18e-3;
  WM_GSP[5] = 29e-3;

/* CP_GSP(j) = CP_GSP_0(j) + CP_GSP_1(j) * T[I] + CP_GSP_2(j) * T[I]^2 */
  CP_GSP_0[1] = 770.0;
  CP_GSP_1[1] = 6.29e-1;
  CP_GSP_2[1] = -1.91e-4;
      
  CP_GSP_0[2] = -100.0;
  CP_GSP_1[2] = 4.4;
  CP_GSP_2[2] = -1.57e-3;
      
  CP_GSP_0[3] = -100.0;
  CP_GSP_1[3] = 4.4e0;
  CP_GSP_2[3] = -1.57e-3;
      
  CP_GSP_0[4] = 1670.0;
  CP_GSP_1[4] = 6.4e-1;
  CP_GSP_2[4] = 0.0;
      
  CP_GSP_0[5] = 950.0;
  CP_GSP_1[5] = 1.88e-1;
  CP_GSP_2[5] = 0.0;

/* Sigma: characteristic length(dA) */
  D_SIGMA[1] = (3.690 + 3.941)/2.0; // 50%CO + 50%CO2
  D_SIGMA[2] = 5.349;  //Benzene
  D_SIGMA[3] = 5.349;
  D_SIGMA[4] = 2.641;
  D_SIGMA[5] = 3.711;
      
/* Epsilon/k: K */
  D_EPSILON_K[1] = (91.7 + 195.2)/2.0; // 50%CO + 50%CO2
  D_EPSILON_K[2] = 412.3;
  D_EPSILON_K[3] = 412.3;
  D_EPSILON_K[4] = 809.1;
  D_EPSILON_K[5] = 78.6;

  for (I = 1; I <= NGSP; I++) {
    for (SPC = 1; SPC <= (I-1); SPC++) {
      WM_AB = 2000.0 / (1.0/WM_GSP[I] + 1.0/WM_GSP[SPC]); //In g/mol
      SIGMA_AB = (D_SIGMA[I] + D_SIGMA[SPC]) / 2.0;
//     EPSILON_K = (D_EPSILON_K(I) * D_EPSILON_K(j))^0.5;        
//     Omega: diffusion collision integral, dimensionless
//     T = 590;
//     T_PRIMA = T / EPSILON_K;
//     OMEGA = 1.06036/(T_PRIMA^0.15610) + 0.193/exp(0.47635*T_PRIMA) ...
//     + 1.03587/exp(1.52996*T_PRIMA) + 1.76474/exp(3.89411*T_PRIMA);
//     D_AB: diffusion coefficient, m2/s
/* Considering pressure: 1 bar
   Include tortuosity */
      D_AB[I][SPC] = 1.0e-4 * 0.00266 / (1.0 * sqrt(WM_AB) * SIGMA_AB*SIGMA_AB);
//     D_AB(i,j) = 1.0e-4 * 0.00266 * T^(3.0/2.0) / ...
//     (PRES * WM_AB^0.5 * SIGMA_AB^2.0 * OMEGA);
/* T and OMEGA are added later */
    }
  }

  VISCOSITY_GAS = 3.0e-5;
  THERM_CONDUC_GAS = 0.0258;
  DIFUSS_FACTOR = 2.0;
  CTE_GAS = 8.314472;

/*******************
  BOUNDARY CONDITIONS
 *******************/

/* Gas pyrolysis; 1:Other gas, 2:Tar1, 3:Tar2, 4:H2O, 5:Dry air */

  Y_GAS_OUT[1] = 0.0;
  Y_GAS_OUT[2] = 0.0;
  Y_GAS_OUT[3] = 0.0;
  Y_GAS_OUT[4] = 0.0;
  Y_GAS_OUT[5] = 1.0;

  T_OUT = 900.0;

  PRES_OUT = 101325.0;

  WM_GAS_OUT = 0.0;
  for (I = 1; I <= NGSP; I++) {
    WM_GAS_OUT = WM_GAS_OUT + Y_GAS_OUT[I]/WM_GSP[I];
  }
  WM_GAS_OUT = 1.0 / WM_GAS_OUT;
    
  for (I = 1; I <= NGSP; I++) {
    X_GAS_OUT[I] = Y_GAS_OUT[I] * WM_GAS_OUT / WM_GSP[I];
  }
  DENSITY_GAS_OUT = PRES_OUT * WM_GAS_OUT / (CTE_GAS * T_OUT);

  F_FLUX = 0.0;

  ALPHA = 50.0;
  BETA = 0.003;
  S_EMISS = 0.85;

/******************
  INITIAL CONDITIONS
 ******************/

  POROSITY_0 = 0.68;

  Y_GAS_0[1] = 0.0;
  Y_GAS_0[2] = 0.0;
  Y_GAS_0[3] = 0.0;
  Y_GAS_0[4] = 0.0;
  Y_GAS_0[5] = 1.0;
      
  Y_SOLID_0[1] = 1.0;
  Y_SOLID_0[2] = 0.0;
      
  T_0 = 300.0;
      
  P_0 = PRES_OUT;
      
  WM_GAS_0 = 0.0;
  for (I = 1; I <= NGSP; I++) {
    WM_GAS_0 = WM_GAS_0 + Y_GAS_0[I]/WM_GSP[I];
  }
  WM_GAS_0 = 1.0 / WM_GAS_0;
    
  DENSITY_GAS_0 =  P_0 * WM_GAS_0 / (CTE_GAS * T_0);

/* Scalex = Scale factor of the x ODE */
  SCALE_SOLID = 500;
  SCALE_T = 1000;
}	  

#include <windows.h>
double getTime(void) {
  SYSTEMTIME st;
  GetSystemTime(&st);
  return (st.wMinute*60 + st.wSecond + (double)st.wMilliseconds/1000);
}
