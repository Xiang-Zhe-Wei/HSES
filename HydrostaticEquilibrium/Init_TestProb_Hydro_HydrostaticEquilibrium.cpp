#include "GAMER.h"

// ===== Hernquist density (Eq. 2) =====
static inline real RhoHernquist(const real r, const real M, const real a) // gamer/include/Typedef.h
{  
   return M*a / ( 2 * M_PI * r * pow( r+a, 3 ) );
}
// =====================================


// https://hackmd.io/mH2qiL4zRii5Pbz6Tn6ZcA?view -- Pressure equation
// ====== Pressure Hernquist ===========
static inline real PressureHernquist(const real r, const real M, const real a){
   const real G = NEWTON_G;          
   const real ap = r + a;
   real term_log = log(ap / r);
   real term1   = a / ap;
   real term2   = a*a / (2.0*ap*ap);
   real term3   = a*a*a / (3.0*ap*ap*ap);
   real term4   = a*a*a*a / (4.0*ap*ap*ap*ap);
   return (G * M * M) / (2.0*M_PI*pow(a,4.0)) * ( term_log - term1 - term2 - term3 - term4 );
}
// =====================================





// problem-specific global variables
// =======================================================================================
static double Blast_Center[3];      // explosion center
static double Mtot;                 // total mass
static double a;                    // scale length
// =======================================================================================

// problem-specific function prototypes
bool Flag_BlastWave( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );




//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

   if ( !OPT__INIT_RESTRICT )
      Aux_Error( ERROR_INFO, "OPT__INIT_RESTRICT must be enabled !!\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  LoadInputTestProb
// Description :  Read problem-specific runtime parameters from Input__TestProb and store them in HDF5 snapshots (Data_*)
//
// Note        :  1. Invoked by SetParameter() to read parameters
//                2. Invoked by Output_DumpData_Total_HDF5() using the function pointer Output_HDF5_InputTest_Ptr to store parameters
//                3. If there is no problem-specific runtime parameter to load, add at least one parameter
//                   to prevent an empty structure in HDF5_Output_t
//                   --> Example:
//                       LOAD_PARA( load_mode, "TestProb_ID", &TESTPROB_ID, TESTPROB_ID, TESTPROB_ID, TESTPROB_ID );
//
// Parameter   :  load_mode      : Mode for loading parameters
//                                 --> LOAD_READPARA    : Read parameters from Input__TestProb
//                                     LOAD_HDF5_OUTPUT : Store parameters in HDF5 snapshots
//                ReadPara       : Data structure for reading parameters (used with LOAD_READPARA)
//                HDF5_InputTest : Data structure for storing parameters in HDF5 snapshots (used with LOAD_HDF5_OUTPUT)
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LoadInputTestProb( const LoadParaMode_t load_mode, ReadPara_t *ReadPara, HDF5_Output_t *HDF5_InputTest )
{

#  ifndef SUPPORT_HDF5
   if ( load_mode == LOAD_HDF5_OUTPUT )   Aux_Error( ERROR_INFO, "please turn on SUPPORT_HDF5 in the Makefile for load_mode == LOAD_HDF5_OUTPUT !!\n" );
#  endif

   if ( load_mode == LOAD_READPARA     &&  ReadPara       == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_READPARA and ReadPara == NULL !!\n" );
   if ( load_mode == LOAD_HDF5_OUTPUT  &&  HDF5_InputTest == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_HDF5_OUTPUT and HDF5_InputTest == NULL !!\n" );

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// --> LOAD_PARA() is defined in "include/TestProb.h"
// ***************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",      &VARIABLE,               DEFAULT,      MIN,              MAX               );
// ***************************************************************************************************************************
   LOAD_PARA( load_mode, "Blast_Center_X",       &Blast_Center[0],       -1.0,          NoMin_double,     amr->BoxSize[0]   );
   LOAD_PARA( load_mode, "Blast_Center_Y",       &Blast_Center[1],       -1.0,          NoMin_double,     amr->BoxSize[1]   );
   LOAD_PARA( load_mode, "Blast_Center_Z",       &Blast_Center[2],       -1.0,          NoMin_double,     amr->BoxSize[2]   );
   LOAD_PARA( load_mode, "Mtot",                 &Mtot,                  -1.0,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "a",                    &a,                     -1.0,          Eps_double,       NoMax_double      );

} // FUNCITON : LoadInputTestProb



//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );

// (1) load the problem-specific runtime parameters
// (1-1) read parameters from Input__TestProb
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

   LoadInputTestProb( LOAD_READPARA, ReadPara, NULL );

   ReadPara->Read( FileName );

   delete ReadPara;

// set the default explosion center
   for (int d=0; d<3; d++)
      if ( Blast_Center[d] < 0.0 )  Blast_Center[d] = 0.5*amr->BoxSize[d];

// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const double End_T_Default    = 5.0e-3;
   const long   End_Step_Default = __INT_MAX__;
    if (END_T < 0.0 || END_T >= End_T_Default) {
      // free falling time
      // ref : https://en.wikipedia.org/wiki/Free-fall_time
      // N is large enough to include most of the mass, then be sure the box size exceeds 2N*a
      const double N_max = 10.0;
      const double R = N_max * a;
      const double m = Mtot * R*R / SQR(R + a);
      const double rho_avg = 3 * m / (4 * M_PI * CUBE(R));
      const double tff     = sqrt(3.0*M_PI / (32.0 * NEWTON_G * rho_avg));
      END_T = 3 * tff;

      // sound cross time
      // ref : https://crossfield.ku.edu/A391_2020A/lec19a.pdf
      const double r = 6.5 * a;           
      const double rho0 = RhoHernquist(r, Mtot, a);
      
      const double P0 = 0.5 * PressureHernquist(r, Mtot, a);
      // const double P0 = 1.0 * PressureHernquist(r, Mtot, a);
      // const double P0 = 1.5 * PressureHernquist(r, Mtot, a);
      const double cs = sqrt( P0 / rho0 );                 
      const double tsc = r / cs;

      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message(stdout, "Sound-crossing time t_sc = %e\n", tsc);
      Aux_Message(stdout, "Free-fall time t_ff = %e\n", tff);
      Aux_Message(stdout, "t_sc / t_ff = %e\n", tsc / tff);
      PRINT_RESET_PARA(END_T, FORMAT_REAL, " (auto-set to ~1 t_ff)");
      Aux_Message( stdout, "=============================================================================\n" );
    }


   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

// (4) make a note
   if ( MPI_Rank == 0 )
   {
//    assuming EOS_GAMMA (must not invoke any EoS routine here since it has not been initialized)
#     if ( EOS != EOS_GAMMA )
      Aux_Message( stderr, "WARNING : the total explosion energy below assumes EOS_GAMMA !!\n" );
#     endif
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID           = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  explosion center          = (%13.7e, %13.7e, %13.7e)\n", Blast_Center[0], Blast_Center[1],
                                                                                       Blast_Center[2] );
      Aux_Message( stdout, "  END_T                     = %13.7e\n", END_T );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter


//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time, const int lv, double AuxArray[] )
{
   const double r = SQRT( SQR(x-Blast_Center[0]) + SQR(y-Blast_Center[1]) + SQR(z-Blast_Center[2]) );
   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;
   const real rho = RhoHernquist(r, Mtot, a);
   
   Pres = 0.5 * PressureHernquist(r, Mtot, a);
   // Pres = 1.0 * PressureHernquist(r, Mtot, a);
   // Pres = 1.5 * PressureHernquist(r, Mtot, a);
   Dens = rho;

   MomX = 0.0;
   MomY = 0.0;
   MomZ = 0.0;
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );   // assuming EoS requires no passive scalars
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;
#  endif // #ifdef SRHD ... else ...

} // FUNCTION : SetGridIC


//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_User_Hydro_Equilibrium
// Description :  AMR flag criterion: refine center to lv<2
//-------------------------------------------------------------------------------------------------------
bool Flag_User_Hydro_Equilibrium(const int i, const int j, const int k, const int lv, const int PID, const double *Threshold) {
    const double dh = amr->dh[lv];
    const double Pos[3] = {
        amr->patch[0][lv][PID]->EdgeL[0] + (i + 0.5)*dh,
        amr->patch[0][lv][PID]->EdgeL[1] + (j + 0.5)*dh,
        amr->patch[0][lv][PID]->EdgeL[2] + (k + 0.5)*dh
    };
    const double Center[3] = {0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2]};
    const double dr[3]     = {Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2]};
    const double Radius    = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);

    // simulation time related
    const double R_ref = 2.5 * a;
    return (Radius < R_ref) && (lv <= 2);
}


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_HydrostaticEquilibrium
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_HydrostaticEquilibrium()
{
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// validate the compilation flags and runtime parameters
   Validate();

#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();

// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr        = SetGridIC;
   Flag_User_Ptr                 = Flag_User_Hydro_Equilibrium;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr     = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO )

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_HydrostaticEquilibrium
