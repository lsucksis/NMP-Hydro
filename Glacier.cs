/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2021/5/18
 * Time: 10:41
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;

namespace NoahMP
{
	/// <summary>
	/// Description of NoahMP_Glacier.
	/// </summary>
	public static class Glacier
	{
		public static double GRAV = 9.80616;
		//acceleration due to gravity (m/s2)
		public static double SB = 5.67E-08;
		//Stefan-Boltzmann constant (w/m2/k4)
		public static double VKC = 0.40;
		//von Karman constant
		public static double TFRZ = 273.16;
		//freezing/melting point (k)
		public static double HSUB = 2.8440E06;
		//latent heat of sublimation (j/kg)
		public static double HVAP = 2.5104E06;
		//latent heat of vaporization (j/kg)
		public static double HFUS = 0.3336E06;
		//latent heat of fusion (j/kg)
		public static double CWAT = 4.188E06;
		//specific heat capacity of water (j/m3/k)
		public static double CICE = 2.094E06;
		//specific heat capacity of ice (j/m3/k)
		public static double CPAIR = 1004.64;
		//heat capacity dry air at const pres (j/kg/k)
		public static double TKWAT = 0.6;
		//thermal conductivity of water (w/m/k)
		public static double TKICE = 2.2;
		//thermal conductivity of ice (w/m/k)
		public static double TKAIR = 0.023;
		//thermal conductivity of air (w/m/k)
		public static double RAIR = 287.04;
		//gas constant for dry air (j/kg/k)
		public static double RW = 461.269;
		//gas constant for  water vapor (j/kg/k)
		public static double DENH2O = 1000;
		//density of water (kg/m3)
		public static double DENICE = 917;
		//density of ice (kg/m3)
		public static double MAXVAL(double[] data)
		{
			double maxv = -1e20;
			for (int i = 0; i < data.Length; i++) {
				if (maxv < data[i])
					maxv = data[i];
			}
			return maxv;
		}
		/// <summary>
		/// liquid water holding capacity for snowpack (m3/m3) (0.03)
		/// </summary>
		public static double SSI = 0.03;
		/// <summary>
		/// snow surface roughness length (m) (0.002)
		/// </summary>
		public static double Z0SNO = 0.002;
		/// <summary>
		/// new snow mass to fully cover old snow (mm)
		/// </summary>
		public static double SWEMX = 1.00;
		//
		
		public static void NOAHMP_GLACIER(
			int ILOC, int JLOC, double COSZ, int NSNOW, int NSOIL, double DT,  // IN : Time/Space/Model-related
			double SFCTMP, double SFCPRS, double UU, double VV, double Q2, double SOLDN,  // IN : Forcing
			double PRCP, double LWDN, double TBOT, double ZLVL, FortDoubleArray FICEOLD, FortDoubleArray ZSOIL,  // IN : Forcing
			ref double QSNOW, ref double SNEQVO, ref double ALBOLD, ref double CM, ref double CH, ref int ISNOW,  // IN/OUT : 
			ref double SNEQV, FortDoubleArray SMC, FortDoubleArray ZSNSO, ref double SNOWH, FortDoubleArray SNICE, FortDoubleArray SNLIQ,  // IN/OUT :
			ref double TG, FortDoubleArray STC, FortDoubleArray SH2O, ref double TAUSS, ref double QSFC,           // IN/OUT : 
			out double  FSA, out double FSR, out double FIRA, out double FSH, out double FGEV, out double SSOIL,  // OUT : 
			out double TRAD, out double EDIR, out double RUNSRF, out double RUNSUB, out double SAG, out double ALBEDO,  // OUT :
			out double QSNBOT, out double  PONDING, out double PONDING1, out double PONDING2, out double T2M, out double Q2E,  // OUT :
			out double  EMISSI, out double  FPICE, out double   CH2B                                 // OUT ://
////ifdef WRF_HYDRO
                   , ref double  sfcheadrt                                            
////endif
		)
		{

// --------------------------------------------------------------------------------------------------
// Initial code: Guo-Yue Niu, Oct. 2007
// Modified to glacier: Michael Barlage, June 2012
// --------------------------------------------------------------------------------------------------
			//implicit none
// --------------------------------------------------------------------------------------------------
// input
//  INTEGER                        , INTENT(IN)    :: ILOC   //grid index
//  INTEGER                        , INTENT(IN)    :: JLOC   //grid index
//  REAL                           , INTENT(IN)    :: COSZ   //cosine solar zenith angle [0-1]
//  INTEGER                        , INTENT(IN)    :: NSNOW  //maximum no. of snow layers        
//  INTEGER                        , INTENT(IN)    :: NSOIL  //no. of soil layers        
//  REAL                           , INTENT(IN)    :: DT     //time step [sec]
//  REAL                           , INTENT(IN)    :: SFCTMP //surface air temperature [K]
//  REAL                           , INTENT(IN)    :: SFCPRS //pressure (pa)
//  REAL                           , INTENT(IN)    :: UU     //wind speed in eastward dir (m/s)
//  REAL                           , INTENT(IN)    :: VV     //wind speed in northward dir (m/s)
//  REAL                           , INTENT(IN)    :: Q2     //mixing ratio (kg/kg) lowest model layer
//  REAL                           , INTENT(IN)    :: SOLDN  //downward shortwave radiation (w/m2)
//  REAL                           , INTENT(IN)    :: PRCP   //precipitation rate (kg m-2 s-1)
//  REAL                           , INTENT(IN)    :: LWDN   //downward longwave radiation (w/m2)
//  REAL                           , INTENT(IN)    :: TBOT   //bottom condition for soil temp. [K]
//  REAL                           , INTENT(IN)    :: ZLVL   //reference height (m)
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(IN)    :: FICEOLD//ice fraction at last timestep
//  REAL, DIMENSION(       1:NSOIL), INTENT(IN)    :: ZSOIL  //layer-bottom depth from soil surf (m)
//
////ifdef WRF_HYDRO
//  REAL                           , INTENT(INOUT)    :: sfcheadrt
////endif
//
//// input/output : need arbitary intial values
//  REAL                           , INTENT(INOUT) :: QSNOW  //snowfall [mm/s]
//  REAL                           , INTENT(INOUT) :: SNEQVO //snow mass at last time step (mm)
//  REAL                           , INTENT(INOUT) :: ALBOLD //snow albedo at last time step (CLASS type)
//  REAL                           , INTENT(INOUT) :: CM     //momentum drag coefficient
//  REAL                           , INTENT(INOUT) :: CH     //sensible heat exchange coefficient
//
//// prognostic variables
//  INTEGER                        , INTENT(INOUT) :: ISNOW  //actual no. of snow layers [-]
//  REAL                           , INTENT(INOUT) :: SNEQV  //snow water eqv. [mm]
//  REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SMC    //soil moisture (ice + liq.) [m3/m3]
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: ZSNSO  //layer-bottom depth from snow surf [m]
//  REAL                           , INTENT(INOUT) :: SNOWH  //snow height [m]
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE  //snow layer ice [mm]
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ  //snow layer liquid water [mm]
//  REAL                           , INTENT(INOUT) :: TG     //ground temperature (k)
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC    //snow/soil temperature [k]
//  REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O   //liquid soil moisture [m3/m3]
//  REAL                           , INTENT(INOUT) :: TAUSS  //non-dimensional snow age
//  REAL                           , INTENT(INOUT) :: QSFC   //mixing ratio at lowest model layer
//
//// output
//  REAL                           , INTENT(OUT)   :: FSA    //total absorbed solar radiation (w/m2)
//  REAL                           , INTENT(OUT)   :: FSR    //total reflected solar radiation (w/m2)
//  REAL                           , INTENT(OUT)   :: FIRA   //total net LW rad (w/m2)  [+ to atm]
//  REAL                           , INTENT(OUT)   :: FSH    //total sensible heat (w/m2) [+ to atm]
//  REAL                           , INTENT(OUT)   :: FGEV   //ground evap heat (w/m2) [+ to atm]
//  REAL                           , INTENT(OUT)   :: SSOIL  //ground heat flux (w/m2)   [+ to soil]
//  REAL                           , INTENT(OUT)   :: TRAD   //surface radiative temperature (k)
//  REAL                           , INTENT(OUT)   :: EDIR   //soil surface evaporation rate (mm/s]
//  REAL                           , INTENT(OUT)   :: RUNSRF //surface runoff [mm/s] 
//  REAL                           , INTENT(OUT)   :: RUNSUB //baseflow (saturation excess) [mm/s]
//  REAL                           , INTENT(OUT)   :: SAG    //solar rad absorbed by ground (w/m2)
//  REAL                           , INTENT(OUT)   :: ALBEDO //surface albedo [-]
//  REAL                           , INTENT(OUT)   :: QSNBOT //snowmelt [mm/s]
//  REAL                           , INTENT(OUT)   :: PONDING//surface ponding [mm]
//  REAL                           , INTENT(OUT)   :: PONDING1//surface ponding [mm]
//  REAL                           , INTENT(OUT)   :: PONDING2//surface ponding [mm]
//  REAL                           , INTENT(OUT)   :: T2M     //2-m air temperature over bare ground part [k]
//  REAL                           , INTENT(OUT)   :: Q2E
//  REAL                           , INTENT(OUT)   :: EMISSI
//  REAL                           , INTENT(OUT)   :: FPICE
//  REAL                           , INTENT(OUT)   :: CH2B
//
//// local
//  INTEGER                                        :: IZ     //do-loop index
//  INTEGER, DIMENSION(-NSNOW+1:NSOIL)             :: IMELT  //phase change index [1-melt; 2-freeze]
//  REAL                                           :: RHOAIR //density air (kg/m3)
//  REAL, DIMENSION(-NSNOW+1:NSOIL)                :: DZSNSO //snow/soil layer thickness [m]
//  REAL                                           :: THAIR  //potential temperature (k)
//  REAL                                           :: QAIR   //specific humidity (kg/kg) (q2/(1+q2))
//  REAL                                           :: EAIR   //vapor pressure air (pa)
//  REAL, DIMENSION(       1:    2)                :: SOLAD  //incoming direct solar rad (w/m2)
//  REAL, DIMENSION(       1:    2)                :: SOLAI  //incoming diffuse solar rad (w/m2)
//  REAL, DIMENSION(       1:NSOIL)                :: SICE   //soil ice content (m3/m3)
//  REAL, DIMENSION(-NSNOW+1:    0)                :: SNICEV //partial volume ice of snow [m3/m3]
//  REAL, DIMENSION(-NSNOW+1:    0)                :: SNLIQV //partial volume liq of snow [m3/m3]
//  REAL, DIMENSION(-NSNOW+1:    0)                :: EPORE  //effective porosity [m3/m3]
//  REAL                                           :: QDEW   //ground surface dew rate [mm/s]
//  REAL                                           :: QVAP   //ground surface evap. rate [mm/s]
//  REAL                                           :: LATHEA //latent heat [j/kg]
//  REAL                                           :: QMELT  //internal pack melt
//  REAL                                           :: SWDOWN //downward solar [w/m2]
//  REAL                                           :: BEG_WB //beginning water for error check
			double ZBOT = -8.0;
//
//  CHARACTER*256 message

// --------------------------------------------------------------------------------------------------
// re-process atmospheric forcing

			double THAIR = 0;
			double QAIR = 0;
			double EAIR = 0;
			double RHOAIR = 0;
			double[] SOLAD = new double[2];
			double[] SOLAI = new double[2];
			FortDoubleArray DZSNSO = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortIntArray IMELT = new FortIntArray(-NSNOW + 1, NSOIL);
			FortDoubleArray SNICEV = new FortDoubleArray(-NSNOW + 1, 0);
			FortDoubleArray SNLIQV = new FortDoubleArray(-NSNOW + 1, 0);
			FortDoubleArray EPORE = new FortDoubleArray(-NSNOW + 1, 0);
			FortDoubleArray SICE = new FortDoubleArray(1, NSOIL);
			double QMELT = 0;
			double LATHEA = 0;
			double SWDOWN = 0;

			ATM_GLACIER(SFCPRS, SFCTMP, Q2, SOLDN, COSZ, out THAIR,  
				out QAIR, out EAIR, out RHOAIR, SOLAD, SOLAI, out SWDOWN);

			double BEG_WB = SNEQV;

// snow/soil layer thickness (m); interface depth: ZSNSO < 0; layer thickness DZSNSO > 0

			//DO IZ = ISNOW+1, NSOIL
			for (int IZ = ISNOW + 1; IZ <= NSOIL; IZ++) {
				if (IZ == ISNOW + 1)
					DZSNSO[IZ] = -ZSNSO[IZ];
				else
					DZSNSO[IZ] = ZSNSO[IZ - 1] - ZSNSO[IZ];
         
			}

// compute energy budget (momentum  energy fluxes and phase changes) 

			ENERGY_GLACIER(NSNOW, NSOIL, ISNOW, DT, QSNOW, RHOAIR,  //in
				EAIR, SFCPRS, QAIR, SFCTMP, LWDN, UU,  //in
				VV, SOLAD, SOLAI, COSZ, ZLVL,          //in
				TBOT, ZBOT, ZSNSO, DZSNSO,                  //in
				ref TG, STC, ref SNOWH, ref SNEQV, ref SNEQVO, SH2O,  //inout
				SMC, SNICE, SNLIQ, ref ALBOLD, ref CM, ref CH,  //inout
				ref TAUSS, ref QSFC,                                  //inout
				IMELT, SNICEV, SNLIQV, EPORE, out QMELT, out PONDING,  //out
				out SAG, out FSA, out FSR, out FIRA, out FSH, out FGEV,  //out
				out TRAD, out T2M, out SSOIL, out LATHEA, out Q2E, out EMISSI, out CH2B);   //out
			for (int i = 1; i <= NSOIL; i++) {
				SICE[i] = Math.Max(0.0, SMC[i] - SH2O[i]);
			}
			SNEQVO = SNEQV;

			double QVAP = Math.Max(FGEV / LATHEA, 0);       // positive part of fgev [mm/s] > 0
			double QDEW = Math.Abs(Math.Min(FGEV / LATHEA, 0));  // negative part of fgev [mm/s] > 0
			EDIR = QVAP - QDEW;

// compute water budgets (water storages, ET components, and runoff)

			WATER_GLACIER(NSNOW, NSOIL, IMELT, DT, PRCP, SFCTMP,  //in
				QVAP, QDEW, FICEOLD, ZSOIL,                  //in
				ref ISNOW, ref SNOWH, ref SNEQV, SNICE, SNLIQ, STC,  //inout
				DZSNSO, SH2O, SICE, ref PONDING, ZSNSO,          //inout
				out RUNSRF, out RUNSUB, out QSNOW, out PONDING1,          //out
				out PONDING2, out QSNBOT, out FPICE                             //out
//ifdef WRF_HYDRO
                        , ref sfcheadrt                     
//endif
			);

			if (MAXVAL(SICE.data) < 0.0001) {
				Console.WriteLine("GLACIER HAS MELTED AT:" + ILOC + JLOC + " ARE YOU SURE THIS SHOULD BE A GLACIER POINT?");
				//WRITE(message,*) "GLACIER HAS MELTED AT:",ILOC,JLOC," ARE YOU SURE THIS SHOULD BE A GLACIER POINT?"
				//CALL wrf_debug(10,TRIM(message))
			}
     
// water and energy balance check

			ERROR_GLACIER(ILOC, JLOC, SWDOWN, FSA, FSR, FIRA, 
				FSH, FGEV, SSOIL, SAG, PRCP, EDIR, 
				RUNSRF, RUNSUB, SNEQV, DT, BEG_WB);

			if (SNOWH <= 1E-6 || SNEQV <= 1E-3) {
				SNOWH = 0.0;
				SNEQV = 0.0;
			}

			if (SWDOWN != 0)
				ALBEDO = FSR / SWDOWN;
			else
				ALBEDO = -999.9;
   
    

		}
//END SUBROUTINE NOAHMP_GLACIER
		public static void ATM_GLACIER(double SFCPRS, double SFCTMP, double Q2, double SOLDN, double COSZ,
			out double THAIR, out double QAIR, out double	EAIR, out double RHOAIR, double[] SOLAD, double[] SOLAI,
			out double SWDOWN)
		{
// --------------------------------------------------------------------------------------------------
// re-process atmospheric forcing
// --------------------------------------------------------------------------------------------------

// --------------------------------------------------------------------------------------------------
// inputs

//  REAL                          , INTENT(IN)  :: SFCPRS //pressure (pa)
//  REAL                          , INTENT(IN)  :: SFCTMP //surface air temperature [k]
//  REAL                          , INTENT(IN)  :: Q2     //mixing ratio (kg/kg)
//  REAL                          , INTENT(IN)  :: SOLDN  //downward shortwave radiation (w/m2)
//  REAL                          , INTENT(IN)  :: COSZ   //cosine solar zenith angle [0-1]

// outputs

//  REAL                          , INTENT(OUT) :: THAIR  //potential temperature (k)
//  REAL                          , INTENT(OUT) :: QAIR   //specific humidity (kg/kg) (q2/(1+q2))
//  REAL                          , INTENT(OUT) :: EAIR   //vapor pressure air (pa)
//  REAL, DIMENSION(       1:   2), INTENT(OUT) :: SOLAD  //incoming direct solar radiation (w/m2)
//  REAL, DIMENSION(       1:   2), INTENT(OUT) :: SOLAI  //incoming diffuse solar radiation (w/m2)
//  REAL                          , INTENT(OUT) :: RHOAIR //density air (kg/m3)
//  REAL                          , INTENT(OUT) :: SWDOWN //downward solar filtered by sun angle [w/m2]

//locals

//  REAL                                        :: PAIR   //atm bottom level pressure (pa)
// --------------------------------------------------------------------------------------------------

			double PAIR = SFCPRS;                   // atm bottom level pressure (pa)
			THAIR = SFCTMP * Math.Pow(SFCPRS / PAIR, RAIR / CPAIR);
//       QAIR   = Q2 / (1.0+Q2)           // mixing ratio to specific humidity [kg/kg]
			QAIR = Q2;                       // In WRF, driver converts to specific humidity

			EAIR = QAIR * SFCPRS / (0.622 + 0.378 * QAIR);
			RHOAIR = (SFCPRS - 0.378 * EAIR) / (RAIR * SFCTMP);
			if (COSZ <= 0)
				SWDOWN = 0;
			else
				SWDOWN = SOLDN;
       

			SOLAD[0] = SWDOWN * 0.7 * 0.5;     // direct  vis
			SOLAD[1] = SWDOWN * 0.7 * 0.5;     // direct  nir
			SOLAI[0] = SWDOWN * 0.3 * 0.5;     // diffuse vis
			SOLAI[1] = SWDOWN * 0.3 * 0.5;     // diffuse nir

		}
//END SUBROUTINE ATM_GLACIER
		// ==================================================================================================
		// --------------------------------------------------------------------------------------------------
		public static void ENERGY_GLACIER(int NSNOW, int NSOIL, int ISNOW, double DT, double QSNOW, double RHOAIR,  //in
			double EAIR, double SFCPRS, double QAIR, double SFCTMP, double LWDN, double UU,  //in
			double VV, double[] SOLAD, double[] SOLAI, double COSZ, double ZREF,          //in
			double TBOT, double ZBOT, FortDoubleArray ZSNSO, FortDoubleArray DZSNSO,                  //in
			ref double  TG, FortDoubleArray STC, ref double SNOWH, ref double SNEQV, ref double SNEQVO, FortDoubleArray SH2O,  //inout
			FortDoubleArray   SMC, FortDoubleArray SNICE, FortDoubleArray SNLIQ, ref double ALBOLD, ref double CM, ref double CH,  //inout
			ref double  TAUSS, ref double QSFC,                                  //inout
			FortIntArray   IMELT, FortDoubleArray SNICEV, FortDoubleArray SNLIQV, FortDoubleArray EPORE, out double QMELT, out double PONDING,  //out
			out double SAG, out double FSA, out double FSR, out double FIRA, out double FSH, out double FGEV,  //out
			out double TRAD, out double T2M, out double SSOIL, out double LATHEA, out double Q2E, out double EMISSI, out double CH2B)   //out
		{
// --------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------
//  USE NOAHMP_VEG_PARAMETERS
//  USE NOAHMP_RAD_PARAMETERS
// --------------------------------------------------------------------------------------------------
//  IMPLICIT NONE
// --------------------------------------------------------------------------------------------------
// inputs
//  INTEGER                           , INTENT(IN)    :: NSNOW  //maximum no. of snow layers        
//  INTEGER                           , INTENT(IN)    :: NSOIL  //number of soil layers
//  INTEGER                           , INTENT(IN)    :: ISNOW  //actual no. of snow layers
//  REAL                              , INTENT(IN)    :: DT     //time step [sec]
//  REAL                              , INTENT(IN)    :: QSNOW  //snowfall on the ground (mm/s)
//  REAL                              , INTENT(IN)    :: RHOAIR //density air (kg/m3)
//  REAL                              , INTENT(IN)    :: EAIR   //vapor pressure air (pa)
//  REAL                              , INTENT(IN)    :: SFCPRS //pressure (pa)
//  REAL                              , INTENT(IN)    :: QAIR   //specific humidity (kg/kg)
//  REAL                              , INTENT(IN)    :: SFCTMP //air temperature (k)
//  REAL                              , INTENT(IN)    :: LWDN   //downward longwave radiation (w/m2)
//  REAL                              , INTENT(IN)    :: UU     //wind speed in e-w dir (m/s)
//  REAL                              , INTENT(IN)    :: VV     //wind speed in n-s dir (m/s)
//  REAL   , DIMENSION(       1:    2), INTENT(IN)    :: SOLAD  //incoming direct solar rad. (w/m2)
//  REAL   , DIMENSION(       1:    2), INTENT(IN)    :: SOLAI  //incoming diffuse solar rad. (w/m2)
//  REAL                              , INTENT(IN)    :: COSZ   //cosine solar zenith angle (0-1)
//  REAL                              , INTENT(IN)    :: ZREF   //reference height (m)
//  REAL                              , INTENT(IN)    :: TBOT   //bottom condition for soil temp. (k) 
//  REAL                              , INTENT(IN)    :: ZBOT   //depth for TBOT [m]
//  REAL   , DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)    :: ZSNSO  //layer-bottom depth from snow surf [m]
//  REAL   , DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)    :: DZSNSO //depth of snow  soil layer-bottom [m]
//
//// input  output
//  REAL                              , INTENT(INOUT) :: TG     //ground temperature (k)
//  REAL   , DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC    //snow/soil temperature [k]
//  REAL                              , INTENT(INOUT) :: SNOWH  //snow height [m]
//  REAL                              , INTENT(INOUT) :: SNEQV  //snow mass (mm)
//  REAL                              , INTENT(INOUT) :: SNEQVO //snow mass at last time step (mm)
//  REAL   , DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O   //liquid soil moisture [m3/m3]
//  REAL   , DIMENSION(       1:NSOIL), INTENT(INOUT) :: SMC    //soil moisture (ice + liq.) [m3/m3]
//  REAL   , DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE  //snow ice mass (kg/m2)
//  REAL   , DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ  //snow liq mass (kg/m2)
//  REAL                              , INTENT(INOUT) :: ALBOLD //snow albedo at last time step(CLASS type)
//  REAL                              , INTENT(INOUT) :: CM     //momentum drag coefficient
//  REAL                              , INTENT(INOUT) :: CH     //sensible heat exchange coefficient
//  REAL                              , INTENT(INOUT) :: TAUSS  //snow aging factor
//  REAL                              , INTENT(INOUT) :: QSFC   //mixing ratio at lowest model layer
//
//// outputs
//  INTEGER, DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT)   :: IMELT  //phase change index [1-melt; 2-freeze]
//  REAL   , DIMENSION(-NSNOW+1:    0), INTENT(OUT)   :: SNICEV //partial volume ice [m3/m3]
//  REAL   , DIMENSION(-NSNOW+1:    0), INTENT(OUT)   :: SNLIQV //partial volume liq. water [m3/m3]
//  REAL   , DIMENSION(-NSNOW+1:    0), INTENT(OUT)   :: EPORE  //effective porosity [m3/m3]
//  REAL                              , INTENT(OUT)   :: QMELT  //snowmelt [mm/s]
//  REAL                              , INTENT(OUT)   :: PONDING//pounding at ground [mm]
//  REAL                              , INTENT(OUT)   :: SAG    //solar rad. absorbed by ground (w/m2)
//  REAL                              , INTENT(OUT)   :: FSA    //tot. absorbed solar radiation (w/m2)
//  REAL                              , INTENT(OUT)   :: FSR    //tot. reflected solar radiation (w/m2)
//  REAL                              , INTENT(OUT)   :: FIRA   //total net LW. rad (w/m2)   [+ to atm]
//  REAL                              , INTENT(OUT)   :: FSH    //total sensible heat (w/m2) [+ to atm]
//  REAL                              , INTENT(OUT)   :: FGEV   //ground evaporation (w/m2)  [+ to atm]
//  REAL                              , INTENT(OUT)   :: TRAD   //radiative temperature (k)
//  REAL                              , INTENT(OUT)   :: T2M    //2 m height air temperature (k)
//  REAL                              , INTENT(OUT)   :: SSOIL  //ground heat flux (w/m2)   [+ to soil]
//  REAL                              , INTENT(OUT)   :: LATHEA //latent heat vap./sublimation (j/kg)
//  REAL                              , INTENT(OUT)   :: Q2E
//  REAL                              , INTENT(OUT)   :: EMISSI
//  REAL                              , INTENT(OUT)   :: CH2B   //sensible heat conductance, canopy air to ZLVL air (m/s)
//
//
//// local
//  REAL                                              :: UR     //wind speed at height ZLVL (m/s)
//  REAL                                              :: ZLVL   //reference height (m)
//  REAL                                              :: RSURF  //ground surface resistance (s/m)
//  REAL                                              :: ZPD    //zero plane displacement (m)
//  REAL                                              :: Z0MG   //z0 momentum, ground (m)
//  REAL                                              :: EMG    //ground emissivity
//  REAL                                              :: FIRE   //emitted IR (w/m2)
//  REAL, DIMENSION(-NSNOW+1:NSOIL)                   :: FACT   //temporary used in phase change
//  REAL, DIMENSION(-NSNOW+1:NSOIL)                   :: DF     //thermal conductivity [w/m/k]
//  REAL, DIMENSION(-NSNOW+1:NSOIL)                   :: HCPCT  //heat capacity [j/m3/k]
//  REAL                                              :: GAMMA  //psychrometric constant (pa/k)
//  REAL                                              :: RHSUR  //raltive humidity in surface soil/snow air space (-)

// ---------------------------------------------------------------------------------------------------

// wind speed at reference height: ur >= 1

			double UR = Math.Max(Math.Sqrt(UU * UU + VV * VV), 1);

// roughness length and displacement height

			double Z0MG = NoahMP.Z0SNO;
			double ZPD = SNOWH;

			double ZLVL = ZPD + ZREF;

// Thermal properties of soil, snow, lake, and frozen soil
			FortDoubleArray DF = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortDoubleArray FACT = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortDoubleArray HCPCT = new FortDoubleArray(-NSNOW + 1, NSOIL);
			THERMOPROP_GLACIER(NSOIL, NSNOW, ISNOW, DZSNSO,           //in
				DT, SNOWH, SNICE, SNLIQ,           //in
				DF, HCPCT, SNICEV, SNLIQV, EPORE,  //out
				FACT);                                      //out

// Solar radiation: absorbed  reflected by the ground

			RADIATION_GLACIER(DT, TG, SNEQVO, SNEQV, COSZ,  //in
				QSNOW, SOLAD, SOLAI,                    //in
				ref ALBOLD, ref TAUSS,                             //inout
				out SAG, out FSR, out FSA);                          //out

			// vegetation and ground emissivity
			double EMG = 0.98;
			// soil surface resistance for ground evap.
			double	RHSUR = 1.0;
			double	RSURF = 1.0;

// set psychrometric constant

			LATHEA = HSUB;
			double	GAMMA = CPAIR * SFCPRS / (0.622 * LATHEA);

// Surface temperatures of the ground and energy fluxes

			GLACIER_FLUX(NSOIL, NSNOW, EMG, ISNOW, DF, DZSNSO, Z0MG,  //in
				ZLVL, ZPD, QAIR, SFCTMP, RHOAIR, SFCPRS,  //in
				UR, GAMMA, RSURF, LWDN, RHSUR, SMC,  //in
				EAIR, STC, SAG, SNOWH, LATHEA, SH2O,  //in
				ref CM, ref CH, ref TG, ref QSFC,           //inout
				out FIRA, out FSH, out FGEV, out SSOIL,           //out
				out T2M, out Q2E, out CH2B);                         //out

//energy balance at surface: SAG=(IRB+SHB+EVB+GHB)

			double FIRE = LWDN + FIRA;

			if (FIRE <= 0) {
				//call wrf_error_fatal("STOP in Noah-MP: emitted longwave <0")
			}

			// Compute a net emissivity
			EMISSI = EMG;

			// When we"re computing a TRAD, subtract from the emitted IR the
			// reflected portion of the incoming LWDN, so we"re just
			// considering the IR originating in the canopy/ground system.
    
			TRAD = Math.Pow((FIRE - (1 - EMISSI) * LWDN) / (EMISSI * SB), 0.25);

// 3L snow  4L soil temperatures

			TSNOSOI_GLACIER(NSOIL, NSNOW, ISNOW, DT, TBOT,  //in
				SSOIL, SNOWH, ZBOT, ZSNSO, DF,  //in
				HCPCT,                                      //in
				STC);                                       //inout

// adjusting snow surface temperature
			if (NoahMP.OPT_STC == 2)
			if (SNOWH > 0.05 && TG > TFRZ)
				TG = TFRZ;
     

// Energy released or consumed by snow  frozen soil

			PHASECHANGE_GLACIER(NSNOW, NSOIL, ISNOW, DT, FACT,  //in
				DZSNSO,                                      //in
				STC, SNICE, SNLIQ, ref SNEQV, ref SNOWH,  //inout
				SMC, SH2O,                             //inout
				out QMELT, IMELT, out PONDING);                     //out


		}
// END SUBROUTINE ENERGY_GLACIER
		// ==================================================================================================
		public static void THERMOPROP_GLACIER(int NSOIL, int NSNOW, int ISNOW, FortDoubleArray DZSNSO,  //in
			double DT, double SNOWH, FortDoubleArray SNICE, FortDoubleArray SNLIQ,  //in
			FortDoubleArray DF, FortDoubleArray HCPCT, FortDoubleArray SNICEV, FortDoubleArray SNLIQV, FortDoubleArray EPORE,  //out
			FortDoubleArray FACT)                                       //out
		{
// ------------------------------------------------------------------------------------------------- 
// -------------------------------------------------------------------------------------------------
//  IMPLICIT NONE
//// --------------------------------------------------------------------------------------------------
//// inputs
//  INTEGER                        , INTENT(IN)  :: NSOIL   //number of soil layers
//  INTEGER                        , INTENT(IN)  :: NSNOW   //maximum no. of snow layers        
//  INTEGER                        , INTENT(IN)  :: ISNOW   //actual no. of snow layers
//  REAL                           , INTENT(IN)  :: DT      //time step [s]
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(IN)  :: SNICE   //snow ice mass (kg/m2)
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(IN)  :: SNLIQ   //snow liq mass (kg/m2)
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: DZSNSO  //thickness of snow/soil layers [m]
//  REAL                           , INTENT(IN)  :: SNOWH   //snow height [m]
//
//// outputs
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: DF      //thermal conductivity [w/m/k]
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: HCPCT   //heat capacity [j/m3/k]
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: SNICEV  //partial volume of ice [m3/m3]
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: SNLIQV  //partial volume of liquid water [m3/m3]
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: EPORE   //effective porosity [m3/m3]
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: FACT    //computing energy for phase change
//// --------------------------------------------------------------------------------------------------
//// locals
//
//  INTEGER :: IZ, IZ2
//  REAL, DIMENSION(-NSNOW+1:    0)              :: CVSNO   //volumetric specific heat (j/m3/k)
//  REAL, DIMENSION(-NSNOW+1:    0)              :: TKSNO   //snow thermal conductivity (j/m3/k)
//  REAL                                         :: ZMID    //mid-point soil depth
// --------------------------------------------------------------------------------------------------

// compute snow thermal conductivity and heat capacity
			FortDoubleArray TKSNO = new FortDoubleArray(-NSNOW + 1, 0);
			FortDoubleArray CVSNO = new FortDoubleArray(-NSNOW + 1, 0);
			CSNOW_GLACIER(ISNOW, NSNOW, NSOIL, SNICE, SNLIQ, DZSNSO,  //in
				TKSNO, CVSNO, SNICEV, SNLIQV, EPORE);   //out

			//DO IZ = ISNOW+1, 0
			for (int IZ = -ISNOW + 1; IZ <= 0; IZ++) {
				DF[IZ] = TKSNO[IZ];
				HCPCT[IZ] = CVSNO[IZ];
			}

// compute soil thermal properties (using Noah glacial ice approximations)

			//DO  IZ = 1, NSOIL
			for (int IZ = 1; IZ <= NSOIL; IZ++) {
				double ZMID = 0.5 * (DZSNSO[IZ]);
				//DO IZ2 = 1, IZ-1
				for (int IZ2 = 1; IZ <= IZ - 1; IZ++)
					ZMID += DZSNSO[IZ2];
       
				HCPCT[IZ] = 1E6 * (0.8194 + 0.1309 * ZMID);
				DF[IZ] = 0.32333 + (0.10073 * ZMID);
			}
       
// combine a temporary variable used for melting/freezing of snow and frozen soil

			//DO IZ = ISNOW+1,NSOIL
			for (int IZ = ISNOW + 1; IZ <= NSOIL; IZ++)
				FACT[IZ] = DT / (HCPCT[IZ] * DZSNSO[IZ]);
    
// snow/soil interface

			if (ISNOW == 0)
				DF[1] = (DF[1] * DZSNSO[1] + 0.35 * SNOWH) / (SNOWH + DZSNSO[1]);
			else
				DF[1] = (DF[1] * DZSNSO[1] + DF[0] * DZSNSO[0]) / (DZSNSO[0] + DZSNSO[1]);
    


		}
//END SUBROUTINE THERMOPROP_GLACIER
		// ==================================================================================================
		// --------------------------------------------------------------------------------------------------
		public static void CSNOW_GLACIER(int ISNOW, int NSNOW, int NSOIL, FortDoubleArray SNICE, FortDoubleArray SNLIQ, FortDoubleArray DZSNSO,  //in
			FortDoubleArray TKSNO, FortDoubleArray CVSNO, FortDoubleArray SNICEV, FortDoubleArray SNLIQV, FortDoubleArray EPORE)   //out
		{
// --------------------------------------------------------------------------------------------------
// Snow bulk density,volumetric capacity, and thermal conductivity
//---------------------------------------------------------------------------------------------------
//  IMPLICIT NONE
//---------------------------------------------------------------------------------------------------
// inputs

//  INTEGER,                          INTENT(IN) :: ISNOW  //number of snow layers (-)            
//  INTEGER                        ,  INTENT(IN) :: NSNOW  //maximum no. of snow layers        
//  INTEGER                        ,  INTENT(IN) :: NSOIL  //number of soil layers
//  REAL, DIMENSION(-NSNOW+1:    0),  INTENT(IN) :: SNICE  //snow ice mass (kg/m2)
//  REAL, DIMENSION(-NSNOW+1:    0),  INTENT(IN) :: SNLIQ  //snow liq mass (kg/m2) 
//  REAL, DIMENSION(-NSNOW+1:NSOIL),  INTENT(IN) :: DZSNSO //snow/soil layer thickness [m]
//
//// outputs
//
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: CVSNO  //volumetric specific heat (j/m3/k)
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: TKSNO  //thermal conductivity (w/m/k)
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: SNICEV //partial volume of ice [m3/m3]
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: SNLIQV //partial volume of liquid water [m3/m3]
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: EPORE  //effective porosity [m3/m3]
//
//// locals
//
//  INTEGER :: IZ
//  REAL, DIMENSION(-NSNOW+1:    0) :: BDSNOI  //bulk density of snow(kg/m3)

//---------------------------------------------------------------------------------------------------
// thermal capacity of snow

			//DO IZ = ISNOW+1, 0
			for (int IZ = ISNOW + 1; IZ <= 0; IZ++) {
				SNICEV[IZ] = Math.Min(1, SNICE[IZ] / (DZSNSO[IZ] * DENICE));
				EPORE[IZ] = 1 - SNICEV[IZ];
				SNLIQV[IZ] = Math.Min(EPORE[IZ], SNLIQ[IZ] / (DZSNSO[IZ] * DENH2O));
			}
			FortDoubleArray BDSNOI = new FortDoubleArray(-NSNOW + 1, 0);
			for (int IZ = ISNOW + 1; IZ <= 0; IZ++) {
				BDSNOI[IZ] = (SNICE[IZ] + SNLIQ[IZ]) / DZSNSO[IZ];
				CVSNO[IZ] = CICE * SNICEV[IZ] + CWAT * SNLIQV[IZ];
//      CVSNO[IZ] = 0.525E06                          // constant
			}//enddo

// thermal conductivity of snow
			for (int IZ = ISNOW + 1; IZ <= 0; IZ++) {
				TKSNO[IZ] = 3.2217E-6 * Math.Pow(BDSNOI[IZ], 2);           // Stieglitz(yen,1965)
//    TKSNO[IZ] = 2E-2+2.5E-6*BDSNOI[IZ]*BDSNOI[IZ]   // Anderson, 1976
//    TKSNO[IZ] = 0.35                                // constant
//    TKSNO[IZ] = 2.576E-6*BDSNOI[IZ]**2. + 0.074    // Verseghy (1991)
//    TKSNO[IZ] = 2.22*(BDSNOI[IZ]/1000.)**1.88      // Douvill(Yen, 1981)
			}// ENDDO

		}
//END SUBROUTINE CSNOW_GLACIER
		//===================================================================================================
		public static void RADIATION_GLACIER(double DT, double TG, double SNEQVO, double SNEQV, double COSZ,  //in
			double  QSNOW, double[] SOLAD, double[] SOLAI,                    //in
			ref double  ALBOLD, ref double TAUSS,                             //inout
			out double SAG, out double FSR, out double    FSA)                          //out
		{
// --------------------------------------------------------------------------------------------------
//  IMPLICIT NONE
//// --------------------------------------------------------------------------------------------------
//// input
//  REAL, INTENT(IN)                     :: DT     //time step [s]
//  REAL, INTENT(IN)                     :: TG     //ground temperature (k)
//  REAL, INTENT(IN)                     :: SNEQVO //snow mass at last time step(mm)
//  REAL, INTENT(IN)                     :: SNEQV  //snow mass (mm)
//  REAL, INTENT(IN)                     :: COSZ   //cosine solar zenith angle (0-1)
//  REAL, INTENT(IN)                     :: QSNOW  //snowfall (mm/s)
//  REAL, DIMENSION(1:2)    , INTENT(IN) :: SOLAD  //incoming direct solar radiation (w/m2)
//  REAL, DIMENSION(1:2)    , INTENT(IN) :: SOLAI  //incoming diffuse solar radiation (w/m2)
//
//// inout
//  REAL,                  INTENT(INOUT) :: ALBOLD //snow albedo at last time step (CLASS type)
//  REAL,                  INTENT(INOUT) :: TAUSS  //non-dimensional snow age
//
//// output
//  REAL, INTENT(OUT)                    :: SAG    //solar radiation absorbed by ground (w/m2)
//  REAL, INTENT(OUT)                    :: FSR    //total reflected solar radiation (w/m2)
//  REAL, INTENT(OUT)                    :: FSA    //total absorbed solar radiation (w/m2)
//
//// local
//  INTEGER                              :: IB     //number of radiation bands
//  INTEGER                              :: NBAND  //number of radiation bands
//  REAL                                 :: FAGE   //snow age function (0 - new snow)
//  REAL, DIMENSION(1:2)                 :: ALBSND //snow albedo (direct)
//  REAL, DIMENSION(1:2)                 :: ALBSNI //snow albedo (diffuse)
//  REAL                                 :: ALB    //current CLASS albedo
//  REAL                                 :: Math.Abs    //temporary absorbed rad
//  REAL                                 :: REF    //temporary reflected rad
//  REAL                                 :: FSNO   //snow-cover fraction, = 1 if any snow
//  REAL, DIMENSION(1:2)                 :: ALBICE //albedo land ice: 1=vis, 2=nir
//
//  REAL,PARAMETER :: MPE = 1.E-6

// --------------------------------------------------------------------------------------------------

			int NBAND = 2;
			double[] ALBSND = new double[2];
			double[] ALBSNI = new double[2];
			double[] ALBICE = new double[2];
			ALBICE[0] = 0.80;   //albedo land ice: 1=vis, 2=nir
			ALBICE[1] = 0.55;

// snow age
			double FAGE = 0;
			double ALB = 0;
			SNOW_AGE_GLACIER(DT, TG, SNEQVO, SNEQV, ref TAUSS, out FAGE);

// snow albedos: age even when sun is not present

			if (NoahMP.OPT_ALB == 1)
				SNOWALB_BATS_GLACIER(NBAND, COSZ, FAGE, ALBSND, ALBSNI);
			if (NoahMP.OPT_ALB == 2) {
				SNOWALB_CLASS_GLACIER(NBAND, QSNOW, DT, ref ALB, ALBOLD, ALBSND, ALBSNI);
				ALBOLD = ALB;
			}

// zero summed solar fluxes

			SAG = 0;
			FSA = 0;
			FSR = 0;
   
			double	FSNO = 0.0;
			if (SNEQV > 0.0)
				FSNO = 1.0;

// loop over nband wavebands

			//DO IB = 1, NBAND
			for (int IB = 1; IB <= NBAND; IB++) {
				ALBSND[IB] = ALBICE[IB] * (1 - FSNO) + ALBSND[IB] * FSNO;
				ALBSNI[IB] = ALBICE[IB] * (1 - FSNO) + ALBSNI[IB] * FSNO;
  
// solar radiation absorbed by ground surface

				double ABS = SOLAD[IB] * (1 - ALBSND[IB]) + SOLAI[IB] * (1 - ALBSNI[IB]);
				SAG = SAG + ABS;
				FSA = FSA + ABS;
    
				double REF = SOLAD[IB] * ALBSND[IB] + SOLAI[IB] * ALBSNI[IB];
				FSR = FSR + REF;
    
			}

		}
// END SUBROUTINE RADIATION_GLACIER
		// ==================================================================================================
		public static void SNOW_AGE_GLACIER(double DT, double TG, double SNEQVO, double SNEQV, ref double TAUSS, out double FAGE)
		{
// --------------------------------------------------------------------------------------------------
//  IMPLICIT NONE
// ------------------------ code history ------------------------------------------------------------
// from BATS
// ------------------------ input/output variables --------------------------------------------------
//input
//   REAL, INTENT(IN) :: DT        //main time step (s)
//   REAL, INTENT(IN) :: TG        //ground temperature (k)
//   REAL, INTENT(IN) :: SNEQVO    //snow mass at last time step(mm)
//   REAL, INTENT(IN) :: SNEQV     //snow water per unit ground area (mm)
//
//// inout
//  REAL,  INTENT(INOUT) :: TAUSS  //non-dimensional snow age
//
////output
//   REAL, INTENT(OUT) :: FAGE     //snow age
//
////local
//   REAL            :: TAGE       //total aging effects
//   REAL            :: AGE1       //effects of grain growth due to vapor diffusion
//   REAL            :: AGE2       //effects of grain growth at freezing of melt water
//   REAL            :: AGE3       //effects of soot
//   REAL            :: DELA       //temporary variable
//   REAL            :: SGE        //temporary variable
//   REAL            :: DELS       //temporary variable
//   REAL            :: DELA0      //temporary variable
//   REAL            :: ARG        //temporary variable
// See Yang et al. (1997) J.of Climate for detail.
//---------------------------------------------------------------------------------------------------

			if (SNEQV <= 0.0)
				TAUSS = 0;
			else if (SNEQV > 800)
				TAUSS = 0;
			else {
//          TAUSS = 0.
				double DELA0 = 1E-6 * DT;
				double ARG = 5E3 * (1 / TFRZ - 1 / TG);
				double AGE1 = Math.Exp(ARG);
				double AGE2 = Math.Exp(Math.Min(0, 10 * ARG));
				double AGE3 = 0.3;
				double TAGE = AGE1 + AGE2 + AGE3;
				double DELA = DELA0 * TAGE;
				double DELS = Math.Max(0.0, SNEQV - SNEQVO) / NoahMP.SWEMX;
				double SGE = (TAUSS + DELA) * (1.0 - DELS);
				TAUSS = Math.Max(0, SGE);
			}
			FAGE = TAUSS / (TAUSS + 1);

		}
//END SUBROUTINE SNOW_AGE_GLACIER
		// ==================================================================================================
		// --------------------------------------------------------------------------------------------------
		public static void SNOWALB_BATS_GLACIER(int NBAND, double COSZ, double FAGE, double[] ALBSND, double[] ALBSNI)
		{
// --------------------------------------------------------------------------------------------------
//  IMPLICIT NONE
//// --------------------------------------------------------------------------------------------------
//// input
//
//  INTEGER,INTENT(IN) :: NBAND  //number of waveband classes
//
//  REAL,INTENT(IN) :: COSZ    //cosine solar zenith angle
//  REAL,INTENT(IN) :: FAGE    //snow age correction
//
//// output
//
//  REAL, DIMENSION(1:2),INTENT(OUT) :: ALBSND //snow albedo for direct(1=vis, 2=nir)
//  REAL, DIMENSION(1:2),INTENT(OUT) :: ALBSNI //snow albedo for diffuse
//// ---------------------------------------------------------------------------------------------
//
//  REAL :: FZEN                 //zenith angle correction
//  REAL :: CF1                  //temperary variable
//  REAL :: SL2                  //2.*SL
//  REAL :: SL1                  //1/SL
//  REAL :: SL                   //adjustable parameter
			double C1 = 0.2;  //default in BATS
			double C2 = 0.5;  //default in BATS
//double C1 = 0.2 * 2; // double the default to match Sleepers River"s
//double C2 = 0.5 * 2; // snow surface albedo (double aging effects)
// ---------------------------------------------------------------------------------------------
// zero albedos for all points
			for (int i = 0; i < NBAND; i++) {
				ALBSND[i] = 0;
				ALBSNI[i] = 0;
			}
// when cosz > 0

			double SL = 2.0;
			double SL1 = 1 / SL;
			double SL2 = 2 * SL;
			double CF1 = ((1 + SL1) / (1 + SL2 * COSZ) - SL1);
			double FZEN = Math.Max(CF1, 0);

			ALBSNI[0] = 0.95 * (1 - C1 * FAGE);
			ALBSNI[1] = 0.65 * (1 - C2 * FAGE);

			ALBSND[0] = ALBSNI[0] + 0.4 * FZEN * (1 - ALBSNI[0]);    //  vis direct
			ALBSND[1] = ALBSNI[1] + 0.4 * FZEN * (1 - ALBSNI[1]);    //  nir direct

		}
//END SUBROUTINE SNOWALB_BATS_GLACIER
		// ==================================================================================================
		// --------------------------------------------------------------------------------------------------
		public static void SNOWALB_CLASS_GLACIER(int NBAND, double QSNOW, double DT, ref double ALB, double ALBOLD, double[] ALBSND, double[] ALBSNI)
		{
// --------------------------------------------------------------------------------------------------
//  IMPLICIT NONE
// --------------------------------------------------------------------------------------------------
// input

//  INTEGER,INTENT(IN) :: NBAND  //number of waveband classes
//
//  REAL,INTENT(IN) :: QSNOW     //snowfall (mm/s)
//  REAL,INTENT(IN) :: DT        //time step (sec)
//  REAL,INTENT(IN) :: ALBOLD    //snow albedo at last time step
//
//// in  out
//
//  REAL,                INTENT(INOUT) :: ALB        // 
//// output
//
//  REAL, DIMENSION(1:2),INTENT(OUT) :: ALBSND //snow albedo for direct(1=vis, 2=nir)
//  REAL, DIMENSION(1:2),INTENT(OUT) :: ALBSNI //snow albedo for diffuse
// ---------------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------------
// zero albedos for all points
			for (int i = 0; i < NBAND; i++) {
				ALBSND[i] = 0;
				ALBSNI[i] = 0;
			}
       

// when cosz > 0

			ALB = 0.55 + (ALBOLD - 0.55) * Math.Exp(-0.01 * DT / 3600);

// 1 mm fresh snow(SWE) -- 10mm snow depth, assumed the fresh snow density 100kg/m3
// here assume 1cm snow depth will fully cover the old snow

			if (QSNOW > 0)
				ALB += Math.Min(QSNOW * DT, NoahMP.SWEMX) * (0.84 - ALB) / (NoahMP.SWEMX);         

			ALBSNI[0] = ALB;         // vis diffuse
			ALBSNI[1] = ALB;         // nir diffuse
			ALBSND[0] = ALB;         // vis direct
			ALBSND[1] = ALB;         // nir direct

		}
// END SUBROUTINE SNOWALB_CLASS_GLACIER
		// ==================================================================================================
		public static double TDC(double T, double TFRZ)
		{
			return    Math.Min(50, Math.Max(-50, (T - TFRZ)));
		}
		public static void GLACIER_FLUX(int NSOIL, int NSNOW, double EMG, int ISNOW, FortDoubleArray DF, FortDoubleArray DZSNSO, double Z0M,  //in
			double ZLVL, double ZPD, double QAIR, double SFCTMP, double RHOAIR, double SFCPRS,  //in
			double UR, double GAMMA, double RSURF, double LWDN, double RHSUR, FortDoubleArray SMC,  //in
			double EAIR, FortDoubleArray STC, double SAG, double SNOWH, double LATHEA, FortDoubleArray SH2O,  //in
			ref double CM, ref double CH, ref double TGB, ref double QSFC,           //inout
			out double  IRB, out double SHB, out double EVB, out double GHB,           //out
			out double T2MB, out double Q2B, out double EHB2)                         //out
		{

// --------------------------------------------------------------------------------------------------
// use newton-raphson iteration to solve ground (tg) temperature
// that balances the surface energy budgets for glacier.

// bare soil:
// -SAB + IRB[TG] + SHB[TG] + EVB[TG] + GHB[TG] = 0
// ----------------------------------------------------------------------
//  USE MODULE_MODEL_CONSTANTS
// ----------------------------------------------------------------------
//  IMPLICIT NONE
// ----------------------------------------------------------------------
// input
//  INTEGER, INTENT(IN)                         :: NSNOW  //maximum no. of snow layers        
//  INTEGER, INTENT(IN)                         :: NSOIL  //number of soil layers
//  REAL,                            INTENT(IN) :: EMG    //ground emissivity
//  INTEGER,                         INTENT(IN) :: ISNOW  //actual no. of snow layers
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DF     //thermal conductivity of snow/soil (w/m/k)
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DZSNSO //thickness of snow/soil layers (m)
//  REAL,                            INTENT(IN) :: Z0M    //roughness length, momentum, ground (m)
//  REAL,                            INTENT(IN) :: ZLVL   //reference height (m)
//  REAL,                            INTENT(IN) :: ZPD    //zero plane displacement (m)
//  REAL,                            INTENT(IN) :: QAIR   //specific humidity at height zlvl (kg/kg)
//  REAL,                            INTENT(IN) :: SFCTMP //air temperature at reference height (k)
//  REAL,                            INTENT(IN) :: RHOAIR //density air (kg/m3)
//  REAL,                            INTENT(IN) :: SFCPRS //density air (kg/m3)
//  REAL,                            INTENT(IN) :: UR     //wind speed at height zlvl (m/s)
//  REAL,                            INTENT(IN) :: GAMMA  //psychrometric constant (pa/k)
//  REAL,                            INTENT(IN) :: RSURF  //ground surface resistance (s/m)
//  REAL,                            INTENT(IN) :: LWDN   //atmospheric longwave radiation (w/m2)
//  REAL,                            INTENT(IN) :: RHSUR  //raltive humidity in surface soil/snow air space (-)
//  REAL,                            INTENT(IN) :: EAIR   //vapor pressure air at height (pa)
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: STC    //soil/snow temperature (k)
//  REAL, DIMENSION(       1:NSOIL), INTENT(IN) :: SMC    //soil moisture
//  REAL, DIMENSION(       1:NSOIL), INTENT(IN) :: SH2O   //soil liquid water
//  REAL,                            INTENT(IN) :: SAG    //solar radiation absorbed by ground (w/m2)
//  REAL,                            INTENT(IN) :: SNOWH  //actual snow depth [m]
//  REAL,                            INTENT(IN) :: LATHEA //latent heat of vaporization/subli (j/kg)
//
//// input/output
//  REAL,                         INTENT(INOUT) :: CM     //momentum drag coefficient
//  REAL,                         INTENT(INOUT) :: CH     //sensible heat exchange coefficient
//  REAL,                         INTENT(INOUT) :: TGB    //ground temperature (k)
//  REAL,                         INTENT(INOUT) :: QSFC   //mixing ratio at lowest model layer
//
//// output
//// -SAB + IRB[TG] + SHB[TG] + EVB[TG] + GHB[TG] = 0
//  REAL,                           INTENT(OUT) :: IRB    //net longwave rad (w/m2)   [+ to atm]
//  REAL,                           INTENT(OUT) :: SHB    //sensible heat flux (w/m2) [+ to atm]
//  REAL,                           INTENT(OUT) :: EVB    //latent heat flux (w/m2)   [+ to atm]
//  REAL,                           INTENT(OUT) :: GHB    //ground heat flux (w/m2)  [+ to soil]
//  REAL,                           INTENT(OUT) :: T2MB   //2 m height air temperature (k)
//  REAL,                           INTENT(OUT) :: Q2B    //bare ground heat conductance
//  REAL,                           INTENT(OUT) :: EHB2   //sensible heat conductance for diagnostics
//
//
//// local variables 
//  INTEGER :: NITERB  //number of iterations for surface temperature
//  REAL    :: MPE     //prevents overflow error if division by zero
//  REAL    :: DTG        //change in tg, last iteration (k)
//  INTEGER :: MOZSGN  //number of times MOZ changes sign
//  REAL    :: MOZOLD     //Monin-Obukhov stability parameter from prior iteration
			double FM2 = 0;          //Monin-Obukhov momentum adjustment at 2m
			double FH2 = 0;          //Monin-Obukhov heat adjustment at 2m
			double CH2 = 0;          //Surface exchange at 2m
//double H          //temporary sensible heat flux (w/m2)
//double FV         //friction velocity (m/s)
//double CIR;        //coefficients for ir as function of ts**4
//double CGH;        //coefficients for st as function of ts
			double CSH = 0;        //coefficients for sh as function of ts
			double CEV = 0;        //coefficients for ev as function of esat[ts]
			double CQ2B = 0;       //
//  INTEGER :: ITER    //iteration index
			double Z0H = 0;       //roughness length, sensible heat, ground (m)
			double MOZ = 0;        //Monin-Obukhov stability parameter
			double FM = 0;         //momentum stability correction, weighted by prior iters
			double FH = 0;        //sen heat stability correction, weighted by prior iters
			double RAMB;       //aerodynamic resistance for momentum (s/m)
			double RAHB = 1e-20;       //aerodynamic resistance for sensible heat (s/m)
			double RAWB;       //aerodynamic resistance for water vapor (s/m)
			double ESTG = 0;       //saturation vapor pressure at tg (pa)
			double DESTG = 0;      //d(es)/dt at tg (pa/K)
			double ESATW = 0;      //es for water
			double ESATI = 0;      //es for ice
			double DSATW = 0;      //d(es)/dt at tg (pa/K) for water
			double DSATI = 0;      //d(es)/dt at tg (pa/K) for ice
//  REAL    :: A          //temporary calculation
//  REAL    :: B          //temporary calculation
//  REAL    :: T, TDC     //Kelvin to degree Celsius with limit -50 to +50
//  REAL, DIMENSION(       1:NSOIL) :: SICE   //soil ice

			SHB = -1e30;
			EVB = -1e30;
			IRB = -1e30;
			GHB = -1e30;

// -----------------------------------------------------------------
// initialization variables that do not depend on stability iteration
// -----------------------------------------------------------------
			int NITERB = 5;
			double MPE = 1E-6;
			double DTG = 0;
			double MOZSGN = 0;
			double MOZOLD = 0;
			double H = 0;
			double FV = 0.1;

			double CIR = EMG * SB;
			double CGH = 2 * DF[ISNOW + 1] / DZSNSO[ISNOW + 1];

// -----------------------------------------------------------------
			//loop3: DO ITER = 1, NITERB  // begin stability iteration
			for (int ITER = 0; ITER < NITERB; ITER++) {

				Z0H = Z0M;

//       For now, only allow SFCDIF1 until others can be fixed

				SFCDIF1_GLACIER(ITER, ZLVL, ZPD, Z0H, Z0M,  //in
					QAIR, SFCTMP, H, RHOAIR, MPE, UR,  //in
					ref MOZ, ref MOZSGN, ref FM, ref FH, ref FM2, ref FH2,  //inout
					out FV, out CM, out CH, out CH2);                       //out

				RAMB = Math.Max(1, 1 / (CM * UR));
				RAHB = Math.Max(1, 1 / (CH * UR));
				RAWB = RAHB;

// es and d(es)/dt evaluated at tg

				double	T = TDC(TGB, TFRZ);
				ESAT(T, out ESATW, out ESATI, out DSATW, out DSATI);
				if (T > 0) {
					ESTG = ESATW;
					DESTG = DSATW;
				} else {
					ESTG = ESATI;
					DESTG = DSATI;
				}

				CSH = RHOAIR * CPAIR / RAHB;
				CEV = RHOAIR * CPAIR / GAMMA / (RSURF + RAWB);

// surface fluxes and dtg

				IRB = CIR * Math.Pow(TGB, 4) - EMG * LWDN;
				SHB = CSH * (TGB - SFCTMP);
				EVB = CEV * (ESTG * RHSUR - EAIR);
				GHB = CGH * (TGB - STC[ISNOW + 1]);

				double B = SAG - IRB - SHB - EVB - GHB;
				double A = 4 * CIR * Math.Pow(TGB, 3) + CSH + CEV * DESTG + CGH;
				DTG = B / A;

				IRB = IRB + 4 * CIR * Math.Pow(TGB, 3) * DTG;
				SHB = SHB + CSH * DTG;
				EVB = EVB + CEV * DESTG * DTG;
				GHB = GHB + CGH * DTG;

// update ground surface temperature
				TGB = TGB + DTG;

// for M-O length
				H = CSH * (TGB - SFCTMP);

				T = TDC(TGB, TFRZ);
				ESAT(T, out ESATW, out ESATI, out DSATW, out DSATI);
				if (T > 0)
					ESTG = ESATW;
				else
					ESTG = ESATI;
        
				QSFC = 0.622 * (ESTG * RHSUR) / (SFCPRS - 0.378 * (ESTG * RHSUR));

			}//     END DO loop3 // end stability iteration
// -----------------------------------------------------------------

// if snow on ground and TG > TFRZ: reset TG = TFRZ. reevaluate ground fluxes.
			FortDoubleArray SICE = new FortDoubleArray(1, NSOIL);
			for (int i = 1; i <= NSOIL; i++) {
				SICE[i] = SMC[i] - SH2O[i];
			}
			if (NoahMP.OPT_STC == 1) {
				if ((MAXVAL(SICE.data) > 0.0 || SNOWH > 0.0) && TGB > TFRZ) {
					TGB = TFRZ;
					IRB = CIR * Math.Pow(TGB, 4) - EMG * LWDN;
					SHB = CSH * (TGB - SFCTMP);
					EVB = CEV * (ESTG * RHSUR - EAIR);         //ESTG reevaluate ?
					GHB = SAG - (IRB + SHB + EVB);
				}
			}
// 2m air temperature
			EHB2 = FV * VKC / (Math.Log((2 + Z0H) / Z0H) - FH2);
			CQ2B = EHB2;
			if (EHB2 < 1E-5) {
				T2MB = TGB;
				Q2B = QSFC;
			} else {
				T2MB = TGB - SHB / (RHOAIR * CPAIR) * 1 / EHB2;
				Q2B = QSFC - EVB / (LATHEA * RHOAIR) * (1 / CQ2B + RSURF);
			}

// update CH 
			CH = 1 / RAHB;

		}
//  END SUBROUTINE GLACIER_FLUX
		//  ==================================================================================================
		public static void ESAT(double T, out double ESW, out double ESI, out double DESW, out double DESI)
		{
//---------------------------------------------------------------------------------------------------
// use polynomials to calculate saturation vapor pressure and derivative with
// respect to temperature: over water when t > 0 c and over ice when t <= 0 c
//  IMPLICIT NONE
//---------------------------------------------------------------------------------------------------
// in

//  REAL, intent(in)  :: T              //temperature

//out

//  REAL, intent(out) :: ESW            //saturation vapor pressure over water (pa)
//  REAL, intent(out) :: ESI            //saturation vapor pressure over ice (pa)
//  REAL, intent(out) :: DESW           //d(esat)/dt over water (pa/K)
//  REAL, intent(out) :: DESI           //d(esat)/dt over ice (pa/K)
//
//// local
//
			double A0, A1, A2, A3, A4, A5, A6;  //coefficients for esat over water
			double B0, B1, B2, B3, B4, B5, B6;  //coefficients for esat over ice
			double C0, C1, C2, C3, C4, C5, C6;  //coefficients for dsat over water
			double D0, D1, D2, D3, D4, D5, D6;  //coefficients for dsat over ice

			A0 = 6.107799961;
			A1 = 4.436518521E-01;  
			A2 = 1.428945805E-02;
			A3 = 2.650648471E-04;  
			A4 = 3.031240396E-06;
			A5 = 2.034080948E-08;  
			A6 = 6.136820929E-11;

			B0 = 6.109177956;
			B1 = 5.034698970E-01;  
			B2 = 1.886013408E-02;
			B3 = 4.176223716E-04;  
			B4 = 5.824720280E-06;
			B5 = 4.838803174E-08;  
			B6 = 1.838826904E-10;

			C0 = 4.438099984E-01;
			C1 = 2.857002636E-02;  
			C2 = 7.938054040E-04;
			C3 = 1.215215065E-05;  
			C4 = 1.036561403E-07;
			C5 = 3.532421810e-10;  
			C6 = -7.090244804E-13;

			D0 = 5.030305237E-01;
			D1 = 3.773255020E-02;  
			D2 = 1.267995369E-03;
			D3 = 2.477563108E-05;  
			D4 = 3.005693132E-07;
			D5 = 2.158542548E-09;  
			D6 = 7.131097725E-12;

			ESW = 100 * (A0 + T * (A1 + T * (A2 + T * (A3 + T * (A4 + T * (A5 + T * A6))))));
			ESI = 100 * (B0 + T * (B1 + T * (B2 + T * (B3 + T * (B4 + T * (B5 + T * B6))))));
			DESW = 100 * (C0 + T * (C1 + T * (C2 + T * (C3 + T * (C4 + T * (C5 + T * C6))))));
			DESI = 100 * (D0 + T * (D1 + T * (D2 + T * (D3 + T * (D4 + T * (D5 + T * D6))))));

		}
// END SUBROUTINE ESAT
		// ==================================================================================================

		public static void SFCDIF1_GLACIER(int ITER, double ZLVL, double ZPD, double Z0H, double Z0M,  //in
			double  QAIR, double SFCTMP, double H, double RHOAIR, double MPE, double UR,  //in
			ref double  MOZ, ref double MOZSGN, ref double FM, ref double FH, ref double FM2, ref double FH2,  //inout
			out double  FV, out double CM, out double CH, out double CH2)                  //out
		{
// -------------------------------------------------------------------------------------------------
// computing surface drag coefficient CM for momentum and CH for heat
// -------------------------------------------------------------------------------------------------
//    IMPLICIT NONE
// -------------------------------------------------------------------------------------------------
// inputs
//    INTEGER,              INTENT(IN) :: ITER   //iteration index
//    REAL,                 INTENT(IN) :: ZLVL   //reference height  (m)
//    REAL,                 INTENT(IN) :: ZPD    //zero plane displacement (m)
//    REAL,                 INTENT(IN) :: Z0H    //roughness length, sensible heat, ground (m)
//    REAL,                 INTENT(IN) :: Z0M    //roughness length, momentum, ground (m)
//    REAL,                 INTENT(IN) :: QAIR   //specific humidity at reference height (kg/kg)
//    REAL,                 INTENT(IN) :: SFCTMP //temperature at reference height (k)
//    REAL,                 INTENT(IN) :: H      //sensible heat flux (w/m2) [+ to atm]
//    REAL,                 INTENT(IN) :: RHOAIR //density air (kg/m**3)
//    REAL,                 INTENT(IN) :: MPE    //prevents overflow error if division by zero
//    REAL,                 INTENT(IN) :: UR     //wind speed (m/s)
//
//// in  out
//    REAL,              INTENT(INOUT) :: MOZ    //Monin-Obukhov stability (z/L)
//    INTEGER,           INTENT(INOUT) :: MOZSGN //number of times moz changes sign
//    REAL,              INTENT(INOUT) :: FM     //momentum stability correction, weighted by prior iters
//    REAL,              INTENT(INOUT) :: FH     //sen heat stability correction, weighted by prior iters
//    REAL,              INTENT(INOUT) :: FM2    //sen heat stability correction, weighted by prior iters
//    REAL,              INTENT(INOUT) :: FH2    //sen heat stability correction, weighted by prior iters
//
//// outputs
//    REAL,                INTENT(OUT) :: FV     //friction velocity (m/s)
//    REAL,                INTENT(OUT) :: CM     //drag coefficient for momentum
//    REAL,                INTENT(OUT) :: CH     //drag coefficient for heat
//    REAL,                INTENT(OUT) :: CH2    //drag coefficient for heat
//
//// locals
//    REAL    :: MOZOLD                   //Monin-Obukhov stability parameter from prior iteration
//    REAL    :: TMPCM                    //temporary calculation for CM
//    REAL    :: TMPCH                    //temporary calculation for CH
			double MOL;                     //Monin-Obukhov length (m)
			double TVIR;                    //temporary virtual temperature (k)
			double TMP1, TMP2, TMP3;          //temporary calculation
			double FMNEW;                   //stability correction factor, momentum, for current moz
			double FHNEW;                   //stability correction factor, sen heat, for current moz
			double MOZ2;                   //2/L
//double TMPCM2 ; ;                 //temporary calculation for CM2
//double TMPCH2 ;                  //temporary calculation for CH2
			double FM2NEW;                  //stability correction factor, momentum, for current moz
			double FH2NEW;                   //stability correction factor, sen heat, for current moz
			double TMP12, TMP22, TMP32;        //temporary calculation
//
			double CMFM, CHFH, CM2FM2, CH2FH2;
			FV = 0;
// -------------------------------------------------------------------------------------------------
// Monin-Obukhov stability parameter moz for next iteration

			double	MOZOLD = MOZ;
  
//    if(ZLVL <= ZPD) THEN
//       write(*,*) "critical glacier problem: ZLVL <= ZPD; model stops", zlvl, zpd
//       call wrf_error_fatal("STOP in Noah-MP glacier")
//    ENDIF

			double TMPCM = Math.Log((ZLVL - ZPD) / Z0M);
			double TMPCH = Math.Log((ZLVL - ZPD) / Z0H);
			double TMPCM2 = Math.Log((2.0 + Z0M) / Z0M);
			double TMPCH2 = Math.Log((2.0 + Z0H) / Z0H);

			if (ITER == 1) {
				FV = 0.0;
				MOZ = 0.0;
				MOL = 0.0;
				MOZ2 = 0.0;
			} else {
				TVIR = (1 + 0.61 * QAIR) * SFCTMP;
				TMP1 = VKC * (GRAV / TVIR) * H / (RHOAIR * CPAIR);
				if (Math.Abs(TMP1) <= MPE)
					TMP1 = MPE;
				MOL = -1 * FV * FV * FV / TMP1;
				MOZ = Math.Min((ZLVL - ZPD) / MOL, 1);
				MOZ2 = Math.Min((2.0 + Z0H) / MOL, 1);
			}

// accumulate number of times moz changes sign.

			if (MOZOLD * MOZ < 0)
				MOZSGN = MOZSGN + 1;
			if (MOZSGN >= 2) {
				MOZ = 0;
				FM = 0;
				FH = 0;
				MOZ2 = 0;
				FM2 = 0;
				FH2 = 0;
			}

// evaluate stability-dependent variables using moz from prior iteration
			if (MOZ < 0) {
				TMP1 = Math.Pow(1 - 16 * MOZ, 0.25);
				TMP2 = Math.Log((1 + TMP1 * TMP1) / 2);
				TMP3 = Math.Log((1 + TMP1) / 2);
				FMNEW = 2 * TMP3 + TMP2 - 2 * Math.Atan(TMP1) + 1.5707963;
				FHNEW = 2 * TMP2;

// 2-meter
				TMP12 = Math.Pow(1 - 16 * MOZ2, 0.25);
				TMP22 = Math.Log((1 + TMP12 * TMP12) / 2);
				TMP32 = Math.Log((1 + TMP12) / 2);
				FM2NEW = 2 * TMP32 + TMP22 - 2 * Math.Atan(TMP12) + 1.5707963;
				FH2NEW = 2 * TMP22;
			} else {
				FMNEW = -5 * MOZ;
				FHNEW = FMNEW;
				FM2NEW = -5 * MOZ2;
				FH2NEW = FM2NEW;
			}

// except for first iteration, weight stability factors for previous
// iteration to help avoid flip-flops from one iteration to the next

			if (ITER == 1) {
				FM = FMNEW;
				FH = FHNEW;
				FM2 = FM2NEW;
				FH2 = FH2NEW;
			} else {
				FM = 0.5 * (FM + FMNEW);
				FH = 0.5 * (FH + FHNEW);
				FM2 = 0.5 * (FM2 + FM2NEW);
				FH2 = 0.5 * (FH2 + FH2NEW);
			}

// exchange coefficients

			FH = Math.Min(FH, 0.9 * TMPCH);
			FM = Math.Min(FM, 0.9 * TMPCM);
			FH2 = Math.Min(FH2, 0.9 * TMPCH2);
			FM2 = Math.Min(FM2, 0.9 * TMPCM2);

			CMFM = TMPCM - FM;
			CHFH = TMPCH - FH;
			CM2FM2 = TMPCM2 - FM2;
			CH2FH2 = TMPCH2 - FH2;
			if (Math.Abs(CMFM) <= MPE)
				CMFM = MPE;
			if (Math.Abs(CHFH) <= MPE)
				CHFH = MPE;
			if (Math.Abs(CM2FM2) <= MPE)
				CM2FM2 = MPE;
			if (Math.Abs(CH2FH2) <= MPE)
				CH2FH2 = MPE;
			CM = VKC * VKC / (CMFM * CMFM);
			CH = VKC * VKC / (CMFM * CHFH);
			CH2 = VKC * VKC / (CM2FM2 * CH2FH2);
        
// friction velocity

			FV = UR * Math.Sqrt(CM);
			CH2 = VKC * FV / CH2FH2;

		}
// END SUBROUTINE SFCDIF1_GLACIER
		// ==================================================================================================
		public static void TSNOSOI_GLACIER(int NSOIL, int NSNOW, int ISNOW, double DT, double TBOT,  //in
			double  SSOIL, double SNOWH, double ZBOT, FortDoubleArray ZSNSO, FortDoubleArray DF,  //in
			FortDoubleArray  HCPCT,                                      //in
			FortDoubleArray  STC)                                       //inout
		{
// --------------------------------------------------------------------------------------------------
// Compute snow (up to 3L) and soil (4L) temperature. Note that snow temperatures
// during melting season may exceed melting point (TFRZ) but later in PHASECHANGE
// subroutine the snow temperatures are reset to TFRZ for melting snow.
// --------------------------------------------------------------------------------------------------
//  IMPLICIT NONE
// --------------------------------------------------------------------------------------------------
//input

//    INTEGER,                         INTENT(IN)  :: NSOIL  //no of soil layers (4)
//    INTEGER,                         INTENT(IN)  :: NSNOW  //maximum no of snow layers [3]
//    INTEGER,                         INTENT(IN)  :: ISNOW  //actual no of snow layers
//
//    REAL,                            INTENT(IN)  :: DT     //time step (s)
//    REAL,                            INTENT(IN)  :: TBOT   //
//    REAL,                            INTENT(IN)  :: SSOIL  //ground heat flux (w/m2)
//    REAL,                            INTENT(IN)  :: SNOWH  //snow depth (m)
//    REAL,                            INTENT(IN)  :: ZBOT   //from soil surface (m)
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: ZSNSO  //layer-bot. depth from snow surf.(m)
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: DF     //thermal conductivity
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: HCPCT  //heat capacity (J/m3/k)
//
////input and output
//
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC
//
////local
//
//    INTEGER                                      :: IZ
//    REAL                                         :: ZBOTSNO   //ZBOT from snow surface
//    REAL, DIMENSION(-NSNOW+1:NSOIL)              :: AI, BI, CI, RHSTS
			double EFLXB = 0; //energy influx from soil bottom (w/m2)
//    REAL, DIMENSION(-NSNOW+1:NSOIL)              :: PHI   //light through water (w/m2)

// ----------------------------------------------------------------------
			FortDoubleArray PHI = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortDoubleArray AI = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortDoubleArray BI = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortDoubleArray CI = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortDoubleArray RHSTS = new FortDoubleArray(-NSNOW + 1, NSOIL);

// prescribe solar penetration into ice/snow
			for (int i = ISNOW + 1; i <= NSOIL; i++) {
				PHI[i] = 0;
			}
// adjust ZBOT from soil surface to ZBOTSNO from snow surface

			double ZBOTSNO = ZBOT - SNOWH;   //from snow surface

// compute ice temperatures

			HRT_GLACIER(NSNOW, NSOIL, ISNOW, ZSNSO, 
				STC, TBOT, ZBOTSNO, DF, 
				HCPCT, SSOIL, PHI,            
				AI, BI, CI, RHSTS, 
				out EFLXB);

			HSTEP_GLACIER(NSNOW, NSOIL, ISNOW, DT, 
				AI, BI, CI, RHSTS, 
				STC);

		}
//  END SUBROUTINE TSNOSOI_GLACIER
		// ==================================================================================================
		// ----------------------------------------------------------------------
		public static void HRT_GLACIER(int NSNOW, int  NSOIL, int ISNOW, FortDoubleArray ZSNSO,  //in
			FortDoubleArray STC, double TBOT, double ZBOT, FortDoubleArray DF,  //in
			FortDoubleArray HCPCT, double SSOIL, FortDoubleArray PHI,             //in
			FortDoubleArray AI, FortDoubleArray BI, FortDoubleArray CI, FortDoubleArray RHSTS,  //out
			out double  BOTFLX)                                    //out
		{
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// calculate the right hand side of the time tendency term of the soil
// thermal diffusion equation.  also to compute ( prepare ) the matrix
// coefficients for the tri-diagonal matrix of the implicit time scheme.
// ----------------------------------------------------------------------
//    IMPLICIT NONE
// ----------------------------------------------------------------------
// input

//    INTEGER,                         INTENT(IN)  :: NSOIL  //no of soil layers (4)
//    INTEGER,                         INTENT(IN)  :: NSNOW  //maximum no of snow layers [3]
//    INTEGER,                         INTENT(IN)  :: ISNOW  //actual no of snow layers
//    REAL,                            INTENT(IN)  :: TBOT   //bottom soil temp. at ZBOT (k)
//    REAL,                            INTENT(IN)  :: ZBOT   //depth of lower boundary condition (m)
//                                                           //from soil surface not snow surface
//    REAL,                            INTENT(IN)  :: SSOIL  //ground heat flux (w/m2)
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: ZSNSO  //depth of layer-bottom of snow/soil (m)
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: STC    //snow/soil temperature (k)
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: DF     //thermal conductivity [w/m/k]
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: HCPCT  //heat capacity [j/m3/k]
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: PHI    //light through water (w/m2)
//
//// output
//
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: RHSTS  //right-hand side of the matrix
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: AI     //left-hand side coefficient
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: BI     //left-hand side coefficient
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: CI     //left-hand side coefficient
			BOTFLX = 0; //energy influx from soil bottom (w/m2)
//
//// local
//
//    INTEGER                                      :: K
//    REAL, DIMENSION(-NSNOW+1:NSOIL)              :: DDZ
//    REAL, DIMENSION(-NSNOW+1:NSOIL)              :: DENOM
//    REAL, DIMENSION(-NSNOW+1:NSOIL)              :: DTSDZ
//    REAL, DIMENSION(-NSNOW+1:NSOIL)              :: EFLUX
			double TEMP1 = 0;
			FortDoubleArray EFLUX = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortDoubleArray DDZ = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortDoubleArray DENOM = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortDoubleArray DTSDZ = new FortDoubleArray(-NSNOW + 1, NSOIL);
// ----------------------------------------------------------------------

			//DO K = ISNOW+1, NSOIL
			for (int K = ISNOW + 1; K <= NSOIL; K++) {
				if (K == ISNOW + 1) {
					DENOM[K] = -ZSNSO[K] * HCPCT[K];
					TEMP1 = -ZSNSO[K + 1];
					DDZ[K] = 2.0 / TEMP1;
					DTSDZ[K] = 2.0 * (STC[K] - STC[K + 1]) / TEMP1;
					EFLUX[K] = DF[K] * DTSDZ[K] - SSOIL - PHI[K];
				} else if (K < NSOIL) {
					DENOM[K] = (ZSNSO[K - 1] - ZSNSO[K]) * HCPCT[K];
					TEMP1 = ZSNSO[K - 1] - ZSNSO[K + 1];
					DDZ[K] = 2.0 / TEMP1;
					DTSDZ[K] = 2.0 * (STC[K] - STC[K + 1]) / TEMP1;
					EFLUX[K] = (DF[K] * DTSDZ[K] - DF[K - 1] * DTSDZ[K - 1]) - PHI[K];
				} else if (K == NSOIL) {
					DENOM[K] = (ZSNSO[K - 1] - ZSNSO[K]) * HCPCT[K];
					TEMP1 = ZSNSO[K - 1] - ZSNSO[K];
					if (NoahMP.OPT_TBOT == 1)
						BOTFLX = 0;
           
					if (NoahMP.OPT_TBOT == 2) {
						DTSDZ[K] = (STC[K] - TBOT) / (0.5 * (ZSNSO[K - 1] + ZSNSO[K]) - ZBOT);
						BOTFLX = -DF[K] * DTSDZ[K];
					}
					EFLUX[K] = (-BOTFLX - DF[K - 1] * DTSDZ[K - 1]) - PHI[K];
				}
			}

			//DO K = ISNOW+1, NSOIL
			for (int K = ISNOW + 1; K <= NSOIL; K++) {
				if (K == ISNOW + 1) {
					AI[K] = 0.0;
					CI[K] = -DF[K] * DDZ[K] / DENOM[K];
					if (NoahMP.OPT_STC == 1) {
						BI[K] = -CI[K];
					}
					if (NoahMP.OPT_STC == 2) {
						BI[K] = -CI[K] + DF[K] / (0.5 * ZSNSO[K] * ZSNSO[K] * HCPCT[K]);
					}
				} else if (K < NSOIL) {
					AI[K] = -DF[K - 1] * DDZ[K - 1] / DENOM[K];
					CI[K] = -DF[K] * DDZ[K] / DENOM[K];
					BI[K] = -(AI[K] + CI[K]);
				} else if (K == NSOIL) {
					AI[K] = -DF[K - 1] * DDZ[K - 1] / DENOM[K];
					CI[K] = 0.0;
					BI[K] = -(AI[K] + CI[K]);
				}
				RHSTS[K] = EFLUX[K] / (-DENOM[K]);
			}

		}
//  END SUBROUTINE HRT_GLACIER
		// ==================================================================================================
		// ----------------------------------------------------------------------
		public static void HSTEP_GLACIER(int NSNOW, int NSOIL, int ISNOW, double DT,   //in
			FortDoubleArray  AI, FortDoubleArray BI, FortDoubleArray CI, FortDoubleArray RHSTS,   //inout
			FortDoubleArray  STC)                                     //inout
		{
// ----------------------------------------------------------------------
// CALCULATE/UPDATE THE SOIL TEMPERATURE FIELD.
// ----------------------------------------------------------------------
//    implicit none
// ----------------------------------------------------------------------
// input

//    INTEGER,                         INTENT(IN)    :: NSOIL
//    INTEGER,                         INTENT(IN)    :: NSNOW
//    INTEGER,                         INTENT(IN)    :: ISNOW
//    REAL,                            INTENT(IN)    :: DT
//
//// output  input
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: AI
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: BI
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: CI
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: RHSTS
//
//// local
//    INTEGER                                        :: K
//    REAL, DIMENSION(-NSNOW+1:NSOIL)                :: RHSTSIN
//    REAL, DIMENSION(-NSNOW+1:NSOIL)                :: CIIN
// ----------------------------------------------------------------------
			FortDoubleArray RHSTSIN = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortDoubleArray CIIN = new FortDoubleArray(-NSNOW + 1, NSOIL);
			//DO K = ISNOW+1,NSOIL
			for (int K = ISNOW + 1; K <= NSOIL; K++) {
				RHSTS[K] = RHSTS[K] * DT;
				AI[K] = AI[K] * DT;
				BI[K] = 1 + BI[K] * DT;
				CI[K] = CI[K] * DT;
			}//END DO

// copy values for input variables before call to rosr12

			for (int K = ISNOW + 1; K <= NSOIL; K++) {
				RHSTSIN[K] = RHSTS[K];
				CIIN[K] = CI[K];
			}

// solve the tri-diagonal matrix equation

			ROSR12_GLACIER(CI, AI, BI, CIIN, RHSTSIN, RHSTS, ISNOW + 1, NSOIL, NSNOW);

// update snow  soil temperature

			for (int K = ISNOW + 1; K <= NSOIL; K++) {
				STC[K] = STC[K] + CI[K];
			}

		}
//  END SUBROUTINE HSTEP_GLACIER
		// ==================================================================================================
		public static void ROSR12_GLACIER(FortDoubleArray P, FortDoubleArray A, FortDoubleArray B, FortDoubleArray C, FortDoubleArray D, FortDoubleArray DELTA, int NTOP, int NSOIL, int NSNOW)
		{
// ----------------------------------------------------------------------
// SUBROUTINE ROSR12
// ----------------------------------------------------------------------
// INVERT (SOLVE) THE TRI-DIAGONAL MATRIX PROBLEM SHOWN BELOW:
// //////                                            ////// //////  //////   //////  //////
// //B[1], C[1],  0  ,  0  ,  0  ,   . . .  ,    0   // //      //   //      //
// //A[2], B[2], C[2],  0  ,  0  ,   . . .  ,    0   // //      //   //      //
// // 0  , A[3], B[3], C[3],  0  ,   . . .  ,    0   // //      //   // D[3] //
// // 0  ,  0  , A(4), B(4), C(4),   . . .  ,    0   // // P(4) //   // D(4) //
// // 0  ,  0  ,  0  , A(5), B(5),   . . .  ,    0   // // P(5) //   // D(5) //
// // .                                          .   // //  .   // = //   .  //
// // .                                          .   // //  .   //   //   .  //
// // .                                          .   // //  .   //   //   .  //
// // 0  , . . . , 0 , A(M-2), B(M-2), C(M-2),   0   // //P(M-2)//   //D(M-2)//
// // 0  , . . . , 0 ,   0   , A(M-1), B(M-1), C(M-1)// //P(M-1)//   //D(M-1)//
// // 0  , . . . , 0 ,   0   ,   0   ,  A(M) ,  B(M) // // P(M) //   // D(M) //
// //////                                            ////// //////  //////   //////  //////
// ----------------------------------------------------------------------
//    IMPLICIT NONE
//
//    INTEGER, INTENT(IN)   :: NTOP           
//    INTEGER, INTENT(IN)   :: NSOIL,NSNOW
//    INTEGER               :: K, KK
//
//    REAL, DIMENSION(-NSNOW+1:NSOIL),INTENT(IN):: A, B, D
//    REAL, DIMENSION(-NSNOW+1:NSOIL),INTENT(INOUT):: C,P,DELTA
//FortDoubleArray A=new FortDoubleArray(-NSNOW+1,NSOIL);
//FortDoubleArray B=new FortDoubleArray(-NSNOW+1,NSOIL);
//FortDoubleArray D=new FortDoubleArray(-NSNOW+1,NSOIL);
//FortDoubleArray C=new FortDoubleArray(-NSNOW+1,NSOIL);
//FortDoubleArray P=new FortDoubleArray(-NSNOW+1,NSOIL);
//FortDoubleArray DELTA=new FortDoubleArray(-NSNOW+1,NSOIL);
// ----------------------------------------------------------------------
// INITIALIZE EQN COEF C FOR THE LOWEST SOIL LAYER
// ----------------------------------------------------------------------
			C[NSOIL] = 0.0;
			P[NTOP] = -C[NTOP] / B[NTOP];
// ----------------------------------------------------------------------
// SOLVE THE COEFS FOR THE 1ST SOIL LAYER
// ----------------------------------------------------------------------
			DELTA[NTOP] = D[NTOP] / B[NTOP];
// ----------------------------------------------------------------------
// SOLVE THE COEFS FOR SOIL LAYERS 2 THRU NSOIL
// ----------------------------------------------------------------------
			//DO K = NTOP+1,NSOIL
			for (int K = NTOP + 1; K <= NSOIL; K++) {
				P[K] = -C[K] * (1.0 / (B[K] + A[K] * P[K - 1]));
				DELTA[K] = (D[K] - A[K] * DELTA[K - 1]) * (1.0 / (B[K] + A[K]
				* P[K - 1]));
			}//END DO
// ----------------------------------------------------------------------
// SET P TO DELTA FOR LOWEST SOIL LAYER
// ----------------------------------------------------------------------
			P[NSOIL] = DELTA[NSOIL];
// ----------------------------------------------------------------------
// ADJUST P FOR SOIL LAYERS 2 THRU NSOIL
// ----------------------------------------------------------------------
			for (int K = NTOP + 1; K <= NSOIL; K++) {
				int KK = NSOIL - K + (NTOP - 1) + 1;
				P[KK] = P[KK] * P[KK + 1] + DELTA[KK];
			}
// ----------------------------------------------------------------------
		}
//  END SUBROUTINE ROSR12_GLACIER
		// ----------------------------------------------------------------------
		// ==================================================================================================
		public static void PHASECHANGE_GLACIER(int NSNOW, int NSOIL, int ISNOW, double DT, FortDoubleArray FACT,  //in
			FortDoubleArray  DZSNSO,                                      //in
			FortDoubleArray  STC, FortDoubleArray  SNICE, FortDoubleArray  SNLIQ, ref double  SNEQV, ref double SNOWH,  //inout
			FortDoubleArray SMC, FortDoubleArray   SH2O,                             //inout
			out double QMELT, FortIntArray  IMELT, out double PONDING)                     //out
		{
// ----------------------------------------------------------------------
// melting/freezing of snow water and soil water
// ----------------------------------------------------------------------
//  IMPLICIT NONE
// ----------------------------------------------------------------------
// inputs

//  INTEGER, INTENT(IN)                             :: NSNOW  //maximum no. of snow layers [=3]
//  INTEGER, INTENT(IN)                             :: NSOIL  //No. of soil layers [=4]
//  INTEGER, INTENT(IN)                             :: ISNOW  //actual no. of snow layers [<=3]
//  REAL, INTENT(IN)                                :: DT     //land model time step (sec)
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)     :: FACT   //temporary
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)     :: DZSNSO //snow/soil layer thickness [m]
//
//// inputs/outputs
//
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT)  :: STC    //snow/soil layer temperature [k]
//  REAL, DIMENSION(-NSNOW+1:0)    , INTENT(INOUT)  :: SNICE  //snow layer ice [mm]
//  REAL, DIMENSION(-NSNOW+1:0)    , INTENT(INOUT)  :: SNLIQ  //snow layer liquid water [mm]
//  REAL, INTENT(INOUT)                             :: SNEQV
//  REAL, INTENT(INOUT)                             :: SNOWH
//  REAL, DIMENSION(       1:NSOIL), INTENT(INOUT)  :: SH2O   //soil liquid water [m3/m3]
//  REAL, DIMENSION(       1:NSOIL), INTENT(INOUT)  :: SMC    //total soil water [m3/m3]
//
//// outputs
//  REAL,                               INTENT(OUT) :: QMELT  //snowmelt rate [mm/s]
//  INTEGER, DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: IMELT  //phase change index
//  REAL,                               INTENT(OUT) :: PONDING//snowmelt when snow has no layer [mm]
//
//// local
//
//  INTEGER                         :: J,K         //do loop index
//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: HM        //energy residual [w/m2]
//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: XM        //melting or freezing water [kg/m2]
//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: WMASS0
//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: WICE0 
//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: WLIQ0 
//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: MICE      //soil/snow ice mass [mm]
//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: MLIQ      //soil/snow liquid water mass [mm]
//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: HEATR     //energy residual or loss after melting/freezing
//  REAL                            :: TEMP1     //temporary variables [kg/m2]
//  REAL                            :: PROPOR
//  REAL                            :: XMF       //total latent heat of phase change

// ----------------------------------------------------------------------
// Initialization
			FortDoubleArray HM = new FortDoubleArray(-NSNOW + 1, NSOIL);        //energy residual [w/m2]
			FortDoubleArray XM = new FortDoubleArray(-NSNOW + 1, NSOIL);       //melting or freezing water [kg/m2]
			FortDoubleArray WMASS0 = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortDoubleArray WICE0 = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortDoubleArray WLIQ0 = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortDoubleArray MICE = new FortDoubleArray(-NSNOW + 1, NSOIL);      //soil/snow ice mass [mm]
			FortDoubleArray MLIQ = new FortDoubleArray(-NSNOW + 1, NSOIL);      //soil/snow liquid water mass [mm]
			FortDoubleArray HEATR = new FortDoubleArray(-NSNOW + 1, NSOIL);     //energy residual or loss after melting/freezing
			QMELT = 0;
			PONDING = 0;
			double XMF = 0;

			//DO J = ISNOW+1,0           // all snow layers
			for (int J = ISNOW + 1; J <= 0; J++) {
				MICE[J] = SNICE[J];
				MLIQ[J] = SNLIQ[J];
			}//END DO

			//DO J = 1, NSOIL            // all soil layers
			for (int J = 1; J <= NSOIL; J++) {
				MLIQ[J] = SH2O[J] * DZSNSO[J] * 1000;
				MICE[J] = (SMC[J] - SH2O[J]) * DZSNSO[J] * 1000;
			}

			//DO J = ISNOW+1,NSOIL       // all layers
			for (int J = ISNOW + 1; J <= NSOIL; J++) {
				IMELT[J] = 0;
				HM[J] = 0;
				XM[J] = 0;
				WICE0[J] = MICE[J];
				WLIQ0[J] = MLIQ[J];
				WMASS0[J] = MICE[J] + MLIQ[J];
			}
    
			for (int J = ISNOW + 1; J <= NSOIL; J++) {
				if (MICE[J] > 0 && STC[J] >= TFRZ)   // melting 
         	IMELT[J] = 1;
         
				if (MLIQ[J] > 0 && STC[J] < TFRZ)  // freezing 
         	IMELT[J] = 2;
         

				// If snow exists, but its thickness is not enough to create a layer
				if (ISNOW == 0 && SNEQV > 0 && J == 1) {
					if (STC[J] >= TFRZ)
						IMELT[J] = 1;
				}
			}

// Calculate the energy surplus and loss for melting and freezing

			for (int J = ISNOW + 1; J <= NSOIL; J++) {
				if (IMELT[J] > 0) {
					HM[J] = (STC[J] - TFRZ) / FACT[J];
					STC[J] = TFRZ;
				}

				if (IMELT[J] == 1 && HM[J] < 0) {
					HM[J] = 0;
					IMELT[J] = 0;
				}
				if (IMELT[J] == 2 && HM[J] > 0) {
					HM[J] = 0;
					IMELT[J] = 0;
				}
				XM[J] = HM[J] * DT / HFUS;
			}

// The rate of melting and freezing for snow without a layer, needs more work.

			if (ISNOW == 0 && SNEQV > 0 && XM[1] > 0) {
				double TEMP1 = SNEQV;
				SNEQV = Math.Max(0, TEMP1 - XM[1]);
				double PROPOR = SNEQV / TEMP1;
				SNOWH = Math.Max(0, PROPOR * SNOWH);
				HEATR[1] = HM[1] - HFUS * (TEMP1 - SNEQV) / DT;
				if (HEATR[1] > 0) {
					XM[1] = HEATR[1] * DT / HFUS;
					HM[1] = HEATR[1];
					IMELT[1] = 1;
				} else {
					XM[1] = 0;
					HM[1] = 0;
					IMELT[1] = 0;
				}
				QMELT = Math.Max(0, (TEMP1 - SNEQV)) / DT;
				XMF = HFUS * QMELT;
				PONDING = TEMP1 - SNEQV;
			}

// The rate of melting and freezing for snow and soil

			//DO J = ISNOW+1,NSOIL
			for (int J = ISNOW + 1; J <= NSOIL; J++) {
				if (IMELT[J] > 0 && Math.Abs(HM[J]) > 0) {

					HEATR[J] = 0;
					if (XM[J] > 0) {
						MICE[J] = Math.Max(0, WICE0[J] - XM[J]);
						HEATR[J] = HM[J] - HFUS * (WICE0[J] - MICE[J]) / DT;
					} else if (XM[J] < 0) {
						MICE[J] = Math.Min(WMASS0[J], WICE0[J] - XM[J]);
						HEATR[J] = HM[J] - HFUS * (WICE0[J] - MICE[J]) / DT;
					}

					MLIQ[J] = Math.Max(0, WMASS0[J] - MICE[J]);

					if (Math.Abs(HEATR[J]) > 0) {
						STC[J] = STC[J] + FACT[J] * HEATR[J];
						if (J <= 0) {                       // snow
							if (MLIQ[J] * MICE[J] > 0)
								STC[J] = TFRZ;
						}
					}

					if (J > 0)
						XMF = XMF + HFUS * (WICE0[J] - MICE[J]) / DT;

					if (J < 1) {
						QMELT = QMELT + Math.Max(0, (WICE0[J] - MICE[J])) / DT;
					}
				} //end if
			}//    ENDDO
			for (int i = -NSNOW + 1; i <= NSOIL; i++) {
				HEATR[i] = 0.0;
				XM[i] = 0.0;
			}

// Deal with residuals in ice/soil

// FIRST REMOVE EXCESS HEAT BY REDUCING TEMPERATURE OF LAYERS

			//if (ANY(STC(1:4) > TFRZ) && ANY(STC(1:4) < TFRZ)) THEN
			if ((STC[1] > TFRZ || STC[2] > TFRZ || STC[3] > TFRZ || STC[4] > TFRZ) && (STC[1] < TFRZ || STC[2] < TFRZ || STC[3] < TFRZ || STC[4] < TFRZ)) {
				//DO J = 1,NSOIL
				for (int J = 1; J <= NSOIL; J++) {
					if (STC[J] > TFRZ) {
						HEATR[J] = (STC[J] - TFRZ) / FACT[J];
						for (int K = 1; K <= NSOIL; K++) {
							if (J != K && STC[K] < TFRZ && HEATR[J] > 0.1) {
								HEATR[K] = (STC[K] - TFRZ) / FACT[K];
								if (Math.Abs(HEATR[K]) > HEATR[J]) {  // LAYER ABSORBS ALL
									HEATR[K] = HEATR[K] + HEATR[J];
									STC[K] = TFRZ + HEATR[K] * FACT[K];
									HEATR[J] = 0.0;
								} else {
									HEATR[J] = HEATR[J] + HEATR[K];
									HEATR[K] = 0.0;
									STC[K] = TFRZ;
								}
							}
						}
						STC[J] = TFRZ + HEATR[J] * FACT[J];
					}
				}
			}

// NOW REMOVE EXCESS COLD BY INCREASING TEMPERATURE OF LAYERS (MAY NOT BE NECESSARY WITH ABOVE LOOP)

			if ((STC[1] > TFRZ || STC[2] > TFRZ || STC[3] > TFRZ || STC[4] > TFRZ) && (STC[1] < TFRZ || STC[2] < TFRZ || STC[3] < TFRZ || STC[4] < TFRZ)) {
				for (int J = 1; J <= NSOIL; J++) {
					if (STC[J] < TFRZ) {
						HEATR[J] = (STC[J] - TFRZ) / FACT[J];
						for (int K = 1; K <= NSOIL; K++) {
							if (J != K && STC[K] > TFRZ && HEATR[J] < -0.1) {
								HEATR[K] = (STC[K] - TFRZ) / FACT[K];
								if (HEATR[K] > Math.Abs(HEATR[J])) {  // LAYER ABSORBS ALL
									HEATR[K] = HEATR[K] + HEATR[J];
									STC[K] = TFRZ + HEATR[K] * FACT[K];
									HEATR[J] = 0.0;
								} else {
									HEATR[J] = HEATR[J] + HEATR[K];
									HEATR[K] = 0.0;
									STC[K] = TFRZ;
								}
							}
						}
						STC[J] = TFRZ + HEATR[J] * FACT[J];
					}
				}
			}

// NOW REMOVE EXCESS HEAT BY MELTING ICE

			//if (ANY(STC(1:4) > TFRZ) && ANY(MICE(1:4) > 0.)) THEN
			if ((STC[1] > TFRZ || STC[2] > TFRZ || STC[3] > TFRZ || STC[4] > TFRZ) && (MICE[1] > TFRZ || MICE[2] > TFRZ || MICE[3] > TFRZ || MICE[4] > TFRZ)) {
				for (int J = 1; J <= NSOIL; J++) {
					if (STC[J] > TFRZ) {
						HEATR[J] = (STC[J] - TFRZ) / FACT[J];
						XM[J] = HEATR[J] * DT / HFUS;
						for (int K = 1; K <= NSOIL; K++) {
							if (J != K && MICE[K] > 0 && XM[J] > 0.1) {
								if (MICE[K] > XM[J]) {  // LAYER ABSORBS ALL
									MICE[K] = MICE[K] - XM[J];
									XMF = XMF + HFUS * XM[J] / DT;
									STC[K] = TFRZ;
									XM[J] = 0.0;
								} else {
									XM[J] = XM[J] - MICE[K];
									XMF = XMF + HFUS * MICE[K] / DT;
									MICE[K] = 0.0;
									STC[K] = TFRZ;
								}
								MLIQ[K] = Math.Max(0, WMASS0[K] - MICE[K]);
							}
						}
						HEATR[J] = XM[J] * HFUS / DT;
						STC[J] = TFRZ + HEATR[J] * FACT[J];
					}
				}
			}

// NOW REMOVE EXCESS COLD BY FREEZING LIQUID OF LAYERS (MAY NOT BE NECESSARY WITH ABOVE LOOP)

			//if (ANY(STC(1:4) < TFRZ) && ANY(MLIQ(1:4) > 0.)) THEN
			if ((STC[1] < TFRZ || STC[2] < TFRZ || STC[3] < TFRZ || STC[4] < TFRZ) && (MLIQ[1] > TFRZ || MLIQ[2] > TFRZ || MLIQ[3] > TFRZ || MLIQ[4] > TFRZ)) {
				for (int J = 1; J <= NSOIL; J++) {
					if (STC[J] < TFRZ) {
						HEATR[J] = (STC[J] - TFRZ) / FACT[J];
						XM[J] = HEATR[J] * DT / HFUS;
						for (int K = 1; K <= NSOIL; K++) {
							if (J != K && MLIQ[K] > 0 && XM[J] < -0.1) {
								if (MLIQ[K] > Math.Abs(XM[J])) {  // LAYER ABSORBS ALL
									MICE[K] = MICE[K] - XM[J];
									XMF = XMF + HFUS * XM[J] / DT;
									STC[K] = TFRZ;
									XM[J] = 0.0;
								} else {
									XM[J] = XM[J] + MLIQ[K];
									XMF = XMF - HFUS * MLIQ[K] / DT;
									MICE[K] = WMASS0[K];
									STC[K] = TFRZ;
								}
								MLIQ[K] = Math.Max(0, WMASS0[K] - MICE[K]);
							}
						}
						HEATR[J] = XM[J] * HFUS / DT;
						STC[J] = TFRZ + HEATR[J] * FACT[J];
					}
				}
			}

			//DO J = ISNOW+1,0             // snow
			for (int J = ISNOW + 1; J <= 0; J++) {
				SNLIQ[J] = MLIQ[J];
				SNICE[J] = MICE[J];
			}

			//DO J = 1, NSOIL              // soil
			for (int J = 1; J <= NSOIL; J++) {
				SH2O[J] = MLIQ[J] / (1000 * DZSNSO[J]);
				SH2O[J] = Math.Max(0.0, Math.Min(1.0, SH2O[J]));
//       SMC[J]  = (MLIQ[J] + MICE[J]) / (1000. * DZSNSO[J])
				SMC[J] = 1.0;
			}
   
		}
//  END SUBROUTINE PHASECHANGE_GLACIER
		// ==================================================================================================
		public static void WATER_GLACIER(int NSNOW, int NSOIL, FortIntArray IMELT, double DT, double PRCP, double SFCTMP,  //in
			double QVAP, double QDEW, FortDoubleArray FICEOLD, FortDoubleArray ZSOIL,                  //in
			ref  int ISNOW, ref double SNOWH, ref double SNEQV, FortDoubleArray SNICE, FortDoubleArray SNLIQ, FortDoubleArray STC,  //inout
			FortDoubleArray  DZSNSO, FortDoubleArray SH2O, FortDoubleArray SICE, ref double PONDING, FortDoubleArray ZSNSO,          //inout
			out double RUNSRF, out double RUNSUB, out double QSNOW, out double PONDING1,		      //out
			out double PONDING2, out double QSNBOT, out double FPICE                             //out
//ifdef WRF_HYDRO
                            , ref double sfcheadrt                                      
//endif
		)  //out
		{
// ----------------------------------------------------------------------  
// Code history:
// Initial code: Guo-Yue Niu, Oct. 2007
// ----------------------------------------------------------------------
//  implicit none
// ----------------------------------------------------------------------
// input
//  INTEGER,                         INTENT(IN)    :: NSNOW   //maximum no. of snow layers
//  INTEGER,                         INTENT(IN)    :: NSOIL   //no. of soil layers
//  INTEGER, DIMENSION(-NSNOW+1:0) , INTENT(IN)    :: IMELT   //melting state index [1-melt; 2-freeze]
//  REAL,                            INTENT(IN)    :: DT      //main time step (s)
//  REAL,                            INTENT(IN)    :: PRCP    //precipitation (mm/s)
//  REAL,                            INTENT(IN)    :: SFCTMP  //surface air temperature [k]
//  REAL,                            INTENT(IN)    :: QVAP    //soil surface evaporation rate[mm/s]
//  REAL,                            INTENT(IN)    :: QDEW    //soil surface dew rate[mm/s]
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(IN)    :: FICEOLD //ice fraction at last timestep
//  REAL, DIMENSION(       1:NSOIL), INTENT(IN)    :: ZSOIL  //layer-bottom depth from soil surf (m)
//
//// input/output
//  INTEGER,                         INTENT(INOUT) :: ISNOW   //actual no. of snow layers
//  REAL,                            INTENT(INOUT) :: SNOWH   //snow height [m]
//  REAL,                            INTENT(INOUT) :: SNEQV   //snow water eqv. [mm]
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE   //snow layer ice [mm]
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ   //snow layer liquid water [mm]
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC     //snow/soil layer temperature [k]
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO  //snow/soil layer thickness [m]
//  REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O    //soil liquid water content [m3/m3]
//  REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SICE    //soil ice content [m3/m3]
//  REAL                           , INTENT(INOUT) :: PONDING //[mm]
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: ZSNSO   //layer-bottom depth from snow surf [m]
//
//// output
//  REAL,                            INTENT(OUT)   :: RUNSRF  //surface runoff [mm/s] 
//  REAL,                            INTENT(OUT)   :: RUNSUB  //baseflow (sturation excess) [mm/s]
//  REAL,                            INTENT(OUT)   :: QSNOW   //snow at ground srf (mm/s) [+]
//  REAL,                            INTENT(OUT)   :: PONDING1
//  REAL,                            INTENT(OUT)   :: PONDING2
//  REAL,                            INTENT(OUT)   :: QSNBOT  //melting water out of snow bottom [mm/s]
//  REAL,                            INTENT(OUT)   :: FPICE   //precipitation frozen fraction
//
//// local
//  REAL                                           :: QRAIN   //rain at ground srf (mm) [+]
//  REAL                                           :: QSEVA   //soil surface evap rate [mm/s]
//  REAL                                           :: QSDEW   //soil surface dew rate [mm/s]
//  REAL                                           :: QSNFRO  //snow surface frost rate[mm/s]
//  REAL                                           :: QSNSUB  //snow surface sublimation rate [mm/s]
//  REAL                                           :: SNOWHIN //snow depth increasing rate (m/s)
//  REAL                                           :: SNOFLOW //glacier flow [mm/s]
//  REAL                                           :: BDFALL  //density of new snow (mm water/m snow)
//  REAL                                           :: REPLACE //replacement water due to sublimation of glacier
			FortDoubleArray SICE_SAVE = new FortDoubleArray(1, NSOIL);  //soil ice content [m3/m3]
			FortDoubleArray SH2O_SAVE = new FortDoubleArray(1, NSOIL);//  REAL, DIMENSION(       1:NSOIL)                :: SH2O_SAVE  //soil liquid water content [m3/m3]
//  INTEGER :: ILEV
//
////ifdef WRF_HYDRO
//  REAL                           , INTENT(INOUT)    :: sfcheadrt
////endif

// ----------------------------------------------------------------------
// initialize

			double SNOFLOW = 0;
			RUNSUB = 0;
			RUNSRF = 0;
			FPICE = 0;
			for (int i = 1; i <= NSOIL; i++) {
				SICE_SAVE[i] = SICE[i];
				SH2O_SAVE[i] = SH2O[i];
			}
// --------------------------------------------------------------------
// partition precipitation into rain and snow (from CANWATER)

// Jordan (1991)

			if (NoahMP.OPT_SNF == 1) {
				if (SFCTMP > TFRZ + 2.5)
					FPICE = 0;
				else if (SFCTMP <= TFRZ + 0.5)
					FPICE = 1.0;
				else if (SFCTMP <= TFRZ + 2)
					FPICE = 1 - (-54.632 + 0.2 * SFCTMP);
				else
					FPICE = 0.6;
         

			}

			if (NoahMP.OPT_SNF == 2) {
				if (SFCTMP >= TFRZ + 2.2)
					FPICE = 0;
				else
					FPICE = 1.0;
//ENDIF
			}// ENDIF

			if (NoahMP.OPT_SNF == 3) {
				if (SFCTMP >= TFRZ)
					FPICE = 0;
				else
					FPICE = 1.0;
       
			}
//     print*, "fpice: ",fpice

// Hedstrom NR and JW Pomeroy (1998), Hydrol. Processes, 12, 1611-1625
// fresh snow density

			double BDFALL = Math.Min(120, 67.92 + 51.25 * Math.Exp((SFCTMP - TFRZ) / 2.59));

			double QRAIN = PRCP * (1 - FPICE);
			QSNOW = PRCP * FPICE;
			double SNOWHIN = QSNOW / BDFALL;
//     print *, "qrain, qsnow",qrain,qsnow,qrain*dt,qsnow*dt

// sublimation, frost, evaporation, and dew

//     QSNSUB = 0.
//     if (SNEQV > 0.) THEN
//       QSNSUB = Math.Min(QVAP, SNEQV/DT)
//     ENDIF
//     QSEVA = QVAP-QSNSUB

//     QSNFRO = 0.
//     if (SNEQV > 0.) THEN
//        QSNFRO = QDEW
//     ENDIF
//     QSDEW = QDEW - QSNFRO

			double QSNSUB = QVAP;  // send total sublimation/frost to SNOWWATER and deal with it there
			double QSNFRO = QDEW;

//     print *, "qvap",qvap,qvap*dt
//     print *, "qsnsub",qsnsub,qsnsub*dt
//     print *, "qseva",qseva,qseva*dt
//     print *, "qsnfro",qsnfro,qsnfro*dt
//     print *, "qdew",qdew,qdew*dt
//     print *, "qsdew",qsdew,qsdew*dt
//print *, "before snowwater", sneqv,snowh,snice,snliq,sh2o,sice
			QSNBOT = 0;
			PONDING1 = 0;
			PONDING2 = 0;
			SNOWWATER_GLACIER(NSNOW, NSOIL, IMELT, DT, SFCTMP,  //in
				SNOWHIN, QSNOW, QSNFRO, QSNSUB, QRAIN,  //in
				FICEOLD, ZSOIL,                          //in
				ref ISNOW, ref SNOWH, ref SNEQV, SNICE, SNLIQ,  //inout
				SH2O, SICE, STC, DZSNSO, ZSNSO,  //inout
				QSNBOT, SNOFLOW, PONDING1, PONDING2);  //out
//print *, "after snowwater", sneqv,snowh,snice,snliq,sh2o,sice
//print *, "ponding", PONDING,PONDING1,PONDING2

			//PONDING: melting water from snow when there is no layer
    
			RUNSRF = (PONDING + PONDING1 + PONDING2) / DT;

			if (ISNOW == 0)
				RUNSRF = RUNSRF + QSNBOT + QRAIN;
			else
				RUNSRF = RUNSRF + QSNBOT;
    

//ifdef WRF_HYDRO
			RUNSRF = RUNSRF + sfcheadrt / DT;  //sfcheadrt units (mm)
//endif
    
			double REPLACE = 0.0;

			//DO ILEV = 1,NSOIL
			for (int ILEV = 1; ILEV <= NSOIL; ILEV++) {
				REPLACE = REPLACE + DZSNSO[ILEV] * (SICE[ILEV] - SICE_SAVE[ILEV] + SH2O[ILEV] - SH2O_SAVE[ILEV]);
			}
			REPLACE = REPLACE * 1000.0 / DT;     // convert to [mm/s]
			for (int i = 1; i <= NSOIL; i++) {
				SICE[i] = Math.Min(1.0, SICE_SAVE[i]);
				SH2O[i] = 1.0 - SICE[i];
			}
//print *, "replace", replace
    
			// use RUNSUB as a water balancer, SNOFLOW is snow that disappears, REPLACE is
			//   water from below that replaces glacier loss

			RUNSUB = SNOFLOW + REPLACE;

		}
//  END SUBROUTINE WATER_GLACIER
		// ==================================================================================================
		// ----------------------------------------------------------------------
		public static void SNOWWATER_GLACIER(int NSNOW, int NSOIL, FortIntArray IMELT, double DT, double SFCTMP,  //in
			double SNOWHIN, double QSNOW, double QSNFRO, double QSNSUB, double QRAIN,  //in
			FortDoubleArray FICEOLD, FortDoubleArray ZSOIL,                          //in
			ref int  ISNOW, ref double SNOWH, ref double SNEQV, FortDoubleArray SNICE, FortDoubleArray SNLIQ,  //inout
			FortDoubleArray   SH2O, FortDoubleArray   SICE, FortDoubleArray  STC, FortDoubleArray  DZSNSO, FortDoubleArray  ZSNSO,  //inout
			double  QSNBOT, double SNOFLOW, double PONDING1, double PONDING2)  //out
		{
// ----------------------------------------------------------------------
//  IMPLICIT NONE
// ----------------------------------------------------------------------
// input
//  INTEGER,                         INTENT(IN)    :: NSNOW  //maximum no. of snow layers
//  INTEGER,                         INTENT(IN)    :: NSOIL  //no. of soil layers
//  INTEGER, DIMENSION(-NSNOW+1:0) , INTENT(IN)    :: IMELT  //melting state index [0-no melt;1-melt]
//  REAL,                            INTENT(IN)    :: DT     //time step (s)
//  REAL,                            INTENT(IN)    :: SFCTMP //surface air temperature [k]
//  REAL,                            INTENT(IN)    :: SNOWHIN//snow depth increasing rate (m/s)
//  REAL,                            INTENT(IN)    :: QSNOW  //snow at ground srf (mm/s) [+]
//  REAL,                            INTENT(IN)    :: QSNFRO //snow surface frost rate[mm/s]
//  REAL,                            INTENT(IN)    :: QSNSUB //snow surface sublimation rate[mm/s]
//  REAL,                            INTENT(IN)    :: QRAIN  //snow surface rain rate[mm/s]
//  REAL, DIMENSION(-NSNOW+1:0)    , INTENT(IN)    :: FICEOLD//ice fraction at last timestep
//  REAL, DIMENSION(       1:NSOIL), INTENT(IN)    :: ZSOIL  //layer-bottom depth from soil surf (m)
//
//// input  output
//  INTEGER,                         INTENT(INOUT) :: ISNOW  //actual no. of snow layers
//  REAL,                            INTENT(INOUT) :: SNOWH  //snow height [m]
//  REAL,                            INTENT(INOUT) :: SNEQV  //snow water eqv. [mm]
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE  //snow layer ice [mm]
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ  //snow layer liquid water [mm]
//  REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O   //soil liquid moisture (m3/m3)
//  REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SICE   //soil ice moisture (m3/m3)
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC    //snow layer temperature [k]
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO //snow/soil layer thickness [m]
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: ZSNSO  //layer-bottom depth from snow surf [m]
//
//// output
//  REAL,                              INTENT(OUT) :: QSNBOT //melting water out of snow bottom [mm/s]
//  REAL,                              INTENT(OUT) :: SNOFLOW//glacier flow [mm]
//  REAL,                              INTENT(OUT) :: PONDING1
//  REAL,                              INTENT(OUT) :: PONDING2

// local
//  INTEGER :: IZ
//  REAL    :: BDSNOW  //bulk density of snow (kg/m3)
// ----------------------------------------------------------------------
			SNOFLOW = 0.0;
			PONDING1 = 0.0;
			PONDING2 = 0.0;

			SNOWFALL_GLACIER(NSOIL, NSNOW, DT, QSNOW, SNOWHIN,  //in
				SFCTMP,                                  //in
				ISNOW, SNOWH, DZSNSO, STC, SNICE,  //inout
				SNLIQ, ref SNEQV);                           //inout

			if (ISNOW < 0) {        //WHEN MORE THAN ONE LAYER
				COMPACT_GLACIER(NSNOW, NSOIL, DT, STC, SNICE,  //in
					SNLIQ, IMELT, FICEOLD,                  //in
					ref ISNOW, DZSNSO);                           //inout

				COMBINE_GLACIER(NSNOW, NSOIL,                          //in
					ref ISNOW, SH2O, STC, SNICE, SNLIQ,  //inout
					DZSNSO, SICE, ref SNOWH, ref SNEQV,          //inout
					ref PONDING1, ref PONDING2);                  //out

				DIVIDE_GLACIER(NSNOW, NSOIL,                          //in
					ref ISNOW, STC, SNICE, SNLIQ, DZSNSO);  //inout
			}

//SET EMPTY SNOW LAYERS TO ZERO

			//DO IZ = -NSNOW+1, ISNOW
			for (int IZ = -NSNOW + 1; IZ <= ISNOW; IZ++) {
				SNICE[IZ] = 0;
				SNLIQ[IZ] = 0;
				STC[IZ] = 0;
				DZSNSO[IZ] = 0;
				ZSNSO[IZ] = 0;
			}//ENDDO

			SNOWH2O_GLACIER(NSNOW, NSOIL, DT, QSNFRO, QSNSUB,  //in 
				QRAIN,                                  //in
				ref ISNOW, DZSNSO, ref SNOWH, ref SNEQV, SNICE,  //inout
				SNLIQ, SH2O, SICE, STC,          //inout
				ref PONDING1, ref PONDING2,          //inout
				out QSNBOT);                                   //out

//to obtain equilibrium state of snow in glacier region
       
			if (SNEQV > 2000) {   // 2000 mm -> maximum water depth
				double BDSNOW = SNICE[0] / DZSNSO[0];
				SNOFLOW = (SNEQV - 2000);
				SNICE[0] = SNICE[0] - SNOFLOW;
				DZSNSO[0] = DZSNSO[0] - SNOFLOW / BDSNOW;
				SNOFLOW = SNOFLOW / DT;
			}

// sum up snow mass for layered snow

			if (ISNOW != 0) {
				SNEQV = 0;
				for (int IZ = ISNOW + 1; IZ <= 0; IZ++) {
					SNEQV = SNEQV + SNICE[IZ] + SNLIQ[IZ];
				}
			}

// Reset ZSNSO and layer thinkness DZSNSO

			for (int IZ = ISNOW + 1; IZ <= 0; IZ++) {
				DZSNSO[IZ] = -DZSNSO[IZ];
			}

			DZSNSO[1] = ZSOIL[1];
			for (int IZ = 2; IZ <= NSOIL; IZ++) {
				DZSNSO[IZ] = (ZSOIL[IZ] - ZSOIL[IZ - 1]);
			}

			ZSNSO[ISNOW + 1] = DZSNSO[ISNOW + 1];
			for (int IZ = ISNOW + 2; IZ <= NSOIL; IZ++) {
				ZSNSO[IZ] = ZSNSO[IZ - 1] + DZSNSO[IZ];
			}

			//DO IZ = ISNOW+1 ,NSOIL
			for (int IZ = ISNOW + 1; IZ <= NSOIL; IZ++) {
				DZSNSO[IZ] = -DZSNSO[IZ];
			}

		}
// END SUBROUTINE SNOWWATER_GLACIER
		// ==================================================================================================
		public static void SNOWFALL_GLACIER(int NSOIL, int NSNOW, double DT, double QSNOW, double SNOWHIN,  //in
			double SFCTMP,                                   //in
			int ISNOW, double SNOWH, FortDoubleArray DZSNSO, FortDoubleArray STC, FortDoubleArray SNICE,  //inout
			FortDoubleArray SNLIQ, ref double SNEQV)                            //inout
		{
// ----------------------------------------------------------------------
// snow depth and density to account for the new snowfall.
// new values of snow depth  density returned.
// ----------------------------------------------------------------------
//    IMPLICIT NONE
// ----------------------------------------------------------------------
// input

//  INTEGER,                            INTENT(IN) :: NSOIL  //no. of soil layers
//  INTEGER,                            INTENT(IN) :: NSNOW  //maximum no. of snow layers
//  REAL,                               INTENT(IN) :: DT     //main time step (s)
//  REAL,                               INTENT(IN) :: QSNOW  //snow at ground srf (mm/s) [+]
//  REAL,                               INTENT(IN) :: SNOWHIN//snow depth increasing rate (m/s)
//  REAL,                               INTENT(IN) :: SFCTMP //surface air temperature [k]
//
//// input and output
//
//  INTEGER,                         INTENT(INOUT) :: ISNOW  //actual no. of snow layers
//  REAL,                            INTENT(INOUT) :: SNOWH  //snow depth [m]
//  REAL,                            INTENT(INOUT) :: SNEQV  //swow water equivalent [m]
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO //thickness of snow/soil layers (m)
//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC    //snow layer temperature [k]
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE  //snow layer ice [mm]
//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ  //snow layer liquid water [mm]
//
//// local
//
//  INTEGER :: NEWNODE            // 0-no new layers, 1-creating new layers
// ----------------------------------------------------------------------
			int NEWNODE = 0;

// shallow snow / no layer

			if (ISNOW == 0 && QSNOW > 0) {
				SNOWH = SNOWH + SNOWHIN * DT;
				SNEQV = SNEQV + QSNOW * DT;
			}

// creating a new layer
 
			if (ISNOW == 0 && QSNOW > 0 && SNOWH >= 0.05) {
				ISNOW = -1;
				NEWNODE = 1;
				DZSNSO[0] = SNOWH;
				SNOWH = 0;
				STC[0] = Math.Min(273.16, SFCTMP);   // temporary setup
				SNICE[0] = SNEQV;
				SNLIQ[0] = 0;
			}

// snow with layers

			if (ISNOW < 0 && NEWNODE == 0 && QSNOW > 0) {
				SNICE[ISNOW + 1] = SNICE[ISNOW + 1] + QSNOW * DT;
				DZSNSO[ISNOW + 1] = DZSNSO[ISNOW + 1] + SNOWHIN * DT;
			}

// ----------------------------------------------------------------------
		}
//  END SUBROUTINE SNOWFALL_GLACIER
		// ==================================================================================================
		// ----------------------------------------------------------------------
		public static void COMPACT_GLACIER(int NSNOW, int NSOIL, double DT, FortDoubleArray STC, FortDoubleArray SNICE,  //in
			FortDoubleArray SNLIQ, FortIntArray IMELT, FortDoubleArray FICEOLD,                 //in
			ref int ISNOW, FortDoubleArray DZSNSO)                          //inout
		{
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
//  IMPLICIT NONE
// ----------------------------------------------------------------------
// input
//   INTEGER,                         INTENT(IN)    :: NSOIL  //no. of soil layers [ =4]
//   INTEGER,                         INTENT(IN)    :: NSNOW  //maximum no. of snow layers [ =3]
//   INTEGER, DIMENSION(-NSNOW+1:0) , INTENT(IN)    :: IMELT  //melting state index [0-no melt;1-melt]
//   REAL,                            INTENT(IN)    :: DT     //time step (sec)
//   REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)    :: STC    //snow layer temperature [k]
//   REAL, DIMENSION(-NSNOW+1:    0), INTENT(IN)    :: SNICE  //snow layer ice [mm]
//   REAL, DIMENSION(-NSNOW+1:    0), INTENT(IN)    :: SNLIQ  //snow layer liquid water [mm]
//   REAL, DIMENSION(-NSNOW+1:    0), INTENT(IN)    :: FICEOLD//ice fraction at last timestep
//
//// input and output
//   INTEGER,                         INTENT(INOUT) :: ISNOW  // actual no. of snow layers
//   REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO // snow layer thickness [m]
//
//// local
			double C2 = 21e-3;  //[m3/kg] // default 21.e-3
			double C3 = 2.5e-6;   //[1/s]
			double C4 = 0.04;    //[1/k]
			double C5 = 2.0;     //
			double DM = 100.0;   //upper Limit on destructive metamorphism compaction [kg/m3]
			double ETA0 = 0.8e+6; //viscosity coefficient [kg-s/m2]
//                                        //according to Anderson, it is between 0.52e6~1.38e6
//   REAL :: BURDEN //pressure of overlying snow [kg/m2]
			double DDZ1 = 0;   //rate of settling of snow pack due to destructive metamorphism.
			double DDZ2 = 0;   //rate of compaction of snow pack due to overburden.
			double DDZ3 = 0;   //rate of compaction of snow pack due to melt [1/s]
//   REAL :: DEXPF  //EXPF=exp(-c4*(273.15-STC)).
//   REAL :: TD     //STC - TFRZ [K]
//   REAL :: PDZDTC //nodal rate of change in fractional-thickness due to compaction [fraction/s]
//   REAL :: VOID   //void (1 - SNICE - SNLIQ)
//   REAL :: WX     //water mass (ice + liquid) [kg/m2]
			double BI = 0;     //partial density of ice [kg/m3]
			FortDoubleArray FICE = new FortDoubleArray(-NSNOW + 1, 0); //   REAL, DIMENSION(-NSNOW+1:0) :: FICE   //fraction of ice at current time step
//
//   INTEGER  :: J

// ----------------------------------------------------------------------
			double BURDEN = 0.0;

			//DO J = ISNOW+1, 0
			for (int J = ISNOW + 1; J <= 0; J++) {
				double WX = SNICE[J] + SNLIQ[J];
				FICE[J] = SNICE[J] / WX;
				double VOID = 1 - (SNICE[J] / DENICE + SNLIQ[J] / DENH2O) / DZSNSO[J];

				// Allow compaction only for non-saturated node and higher ice lens node.
				if (VOID > 0.001 && SNICE[J] > 0.1) {
					BI = SNICE[J] / DZSNSO[J];
					double TD = Math.Max(0, TFRZ - STC[J]);
					double DEXPF = Math.Exp(-C4 * TD);

					// Settling as a result of destructive metamorphism

					DDZ1 = -C3 * DEXPF;

					if (BI > DM)
						DDZ1 = DDZ1 * Math.Exp(-46.0E-3 * (BI - DM));

					// Liquid water term

					if (SNLIQ[J] > 0.01 * DZSNSO[J])
						DDZ1 = DDZ1 * C5;

					// Compaction due to overburden

					DDZ2 = -(BURDEN + 0.5 * WX) * Math.Exp(-0.08 * TD - C2 * BI) / ETA0; // 0.5*WX -> self-burden

					// Compaction occurring during melt

					if (IMELT[J] == 1) {
						DDZ3 = Math.Max(0, (FICEOLD[J] - FICE[J]) / Math.Max(1E-6, FICEOLD[J]));
						DDZ3 = -DDZ3 / DT;           // sometimes too large
					} else {
						DDZ3 = 0;
					}

					// Time rate of fractional change in DZ (units of s-1)

					double PDZDTC = (DDZ1 + DDZ2 + DDZ3) * DT;
					PDZDTC = Math.Max(-0.5, PDZDTC);

					// The change in DZ due to compaction

					DZSNSO[J] = DZSNSO[J] * (1 + PDZDTC);
				}

				// Pressure of overlying snow

				BURDEN = BURDEN + WX;

			}

		}
//  END SUBROUTINE COMPACT_GLACIER
		// ==================================================================================================
		public static void COMBINE_GLACIER(int NSNOW, int NSOIL,                          //in
			ref int  ISNOW, FortDoubleArray SH2O, FortDoubleArray  STC, FortDoubleArray  SNICE, FortDoubleArray  SNLIQ,  //inout
			FortDoubleArray   DZSNSO, FortDoubleArray  SICE, ref double SNOWH, ref double SNEQV,          //inout
			ref double PONDING1, ref double PONDING2)                  //inout
		{
// ----------------------------------------------------------------------
//    IMPLICIT NONE
// ----------------------------------------------------------------------
// input

//    INTEGER, INTENT(IN)     :: NSNOW                        //maximum no. of snow layers
//    INTEGER, INTENT(IN)     :: NSOIL                        //no. of soil layers
//
//// input and output
//
//    INTEGER,                         INTENT(INOUT) :: ISNOW //actual no. of snow layers
//    REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O  //soil liquid moisture (m3/m3)
//    REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SICE  //soil ice moisture (m3/m3)
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC   //snow layer temperature [k]
//    REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE //snow layer ice [mm]
//    REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ //snow layer liquid water [mm]
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO//snow layer depth [m]
//    REAL,                            INTENT(INOUT) :: SNEQV //snow water equivalent [m]
//    REAL,                            INTENT(INOUT) :: SNOWH //snow depth [m]
//    REAL,                            INTENT(INOUT) :: PONDING1
//    REAL,                            INTENT(INOUT) :: PONDING2
//
//// local variables:
//
			// node indices
//    INTEGER :: ISNOW_OLD             // number of top snow layer
			int MSSI = 0;                 // node index
			int NEIBOR = 0;                // adjacent node selected for combination
//    REAL    :: ZWICE                 // total ice mass in snow
//    REAL    :: ZWLIQ                 // total liquid water in snow
//    REAL    :: DZMIN[3]              // minimum of top snow layer
//    DATA DZMIN /0.045, 0.05, 0.2/
			double[] DZMIN = new double[]{ 0.045, 0.05, 0.2 };
//    DATA DZMIN /0.025, 0.025, 0.1/  // MB: change limit
//-----------------------------------------------------------------------

			int ISNOW_OLD = ISNOW;

			//DO J = ISNOW_OLD+1,0
			for (int J = ISNOW_OLD + 1; J <= 0; J++) {
				if (SNICE[J] <= 0.1) {
					if (J != 0) {
						SNLIQ[J + 1] = SNLIQ[J + 1] + SNLIQ[J];
						SNICE[J + 1] = SNICE[J + 1] + SNICE[J];
					} else {
						if (ISNOW_OLD < -1) {
							SNLIQ[J - 1] = SNLIQ[J - 1] + SNLIQ[J];
							SNICE[J - 1] = SNICE[J - 1] + SNICE[J];
						} else {
							PONDING1 = PONDING1 + SNLIQ[J];       // ISNOW WILL GET SET TO ZERO BELOW
							SNEQV = SNICE[J];                    // PONDING WILL GET ADDED TO PONDING FROM
							SNOWH = DZSNSO[J];                   // PHASECHANGE WHICH SHOULD BE ZERO HERE
							SNLIQ[J] = 0.0;                       // BECAUSE THERE IT WAS ONLY CALCULATED
							SNICE[J] = 0.0;                      // FOR THIN SNOW
							DZSNSO[J] = 0.0;
						}
//                SH2O[1] = SH2O[1]+SNLIQ[J]/(DZSNSO[1]*1000.)
//                SICE[1] = SICE[1]+SNICE[J]/(DZSNSO[1]*1000.)
					}

					// shift all elements above this down by one.
					if (J > ISNOW + 1 && ISNOW < -1) {
						//DO I = J, ISNOW+2, -1
						for (int I = J; I >= ISNOW + 2; I--) {
							STC[I] = STC[I - 1];
							SNLIQ[I] = SNLIQ[I - 1];
							SNICE[I] = SNICE[I - 1];
							DZSNSO[I] = DZSNSO[I - 1];
						}
					}
					ISNOW = ISNOW + 1;
				}
			}

// to conserve water in case of too large surface sublimation

			if (SICE[1] < 0) {
				SH2O[1] = SH2O[1] + SICE[1];
				SICE[1] = 0;
			}

			if (ISNOW == 0)
				return;   // MB: get out if no longer multi-layer

			SNEQV = 0;
			SNOWH = 0;
			double ZWICE = 0;
			double ZWLIQ = 0;

			//DO J = ISNOW+1,0
			for (int J = ISNOW + 1; J <= 0; J++) {
				SNEQV = SNEQV + SNICE[J] + SNLIQ[J];
				SNOWH = SNOWH + DZSNSO[J];
				ZWICE = ZWICE + SNICE[J];
				ZWLIQ = ZWLIQ + SNLIQ[J];
			}

// check the snow depth - all snow gone
// the liquid water assumes ponding on soil surface.

//       if (SNOWH < 0.025 && ISNOW < 0 ) THEN // MB: change limit
			if (SNOWH < 0.05 && ISNOW < 0) {
				ISNOW = 0;
				SNEQV = ZWICE;
				PONDING2 = PONDING2 + ZWLIQ;           // LIMIT OF ISNOW < 0 MEANS INPUT PONDING
				if (SNEQV <= 0)
					SNOWH = 0;            // SHOULD BE ZERO; SEE ABOVE
			}

//       if (SNOWH < 0.05 ) THEN
//          ISNOW  = 0
//          SNEQV = ZWICE
//          SH2O[1] = SH2O[1] + ZWLIQ / (DZSNSO[1] * 1000.)
//          if(SNEQV <= 0.) SNOWH = 0.
//       END if

// check the snow depth - snow layers combined

			if (ISNOW < -1) {

				ISNOW_OLD = ISNOW;
				MSSI = 1;

				//DO I = ISNOW_OLD+1,0
				for (int I = ISNOW_OLD + 1; I <= 0; I++) {
					if (DZSNSO[I] < DZMIN[MSSI]) {

						if (I == ISNOW + 1)
							NEIBOR = I + 1;
						else if (I == 0)
							NEIBOR = I - 1;
						else {
							NEIBOR = I + 1;
							if ((DZSNSO[I - 1] + DZSNSO[I]) < (DZSNSO[I + 1] + DZSNSO[I]))
								NEIBOR = I - 1;
						}

						int J, L;
						// Node l and j are combined and stored as node j.
						if (NEIBOR > I) {
							J = NEIBOR;
							L = I;
						} else {
							J = I;
							L = NEIBOR;
						}

						double dzsnso = DZSNSO[J];
						double snliq = SNLIQ[J];
						double snice = SNICE[J];
						double stc = STC[J];
						COMBO_GLACIER(ref dzsnso, ref snliq, ref snice, 
							ref stc, DZSNSO[L], SNLIQ[L], SNICE[L], STC[L]);                
						DZSNSO[J] = dzsnso;
						SNLIQ[J] = snliq;
						SNICE[J] = snice;
						STC[J] = stc;
	
						// Now shift all elements above this down one.
						if (J - 1 > ISNOW + 1) {
							//DO K = J-1, ISNOW+2, -1
							for (int K = J - 1; K >= ISNOW + 2; K--) {
								STC[K] = STC[K - 1];
								SNICE[K] = SNICE[K - 1];
								SNLIQ[K] = SNLIQ[K - 1];
								DZSNSO[K] = DZSNSO[K - 1];
							}
						}

						// Decrease the number of snow layers
						ISNOW = ISNOW + 1;
						if (ISNOW >= -1)
							break;//EXIT
					} else {
						// The layer thickness is greater than the prescribed minimum value
						MSSI = MSSI + 1;

					}
				}

			}
       

		}
//  END SUBROUTINE COMBINE_GLACIER
		// ==================================================================================================

		// ----------------------------------------------------------------------
		public static void COMBO_GLACIER(ref double DZ, ref double WLIQ, ref double WICE, ref double T, double DZ2, double WLIQ2, double WICE2, double T2)
		{
// ----------------------------------------------------------------------
//    IMPLICIT NONE
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------s
// input

//    REAL, INTENT(IN)    :: DZ2   //nodal thickness of 2 elements being combined [m]
//    REAL, INTENT(IN)    :: WLIQ2 //liquid water of element 2 [kg/m2]
//    REAL, INTENT(IN)    :: WICE2 //ice of element 2 [kg/m2]
//    REAL, INTENT(IN)    :: T2    //nodal temperature of element 2 [k]
//    REAL, INTENT(INOUT) :: DZ    //nodal thickness of 1 elements being combined [m]
//    REAL, INTENT(INOUT) :: WLIQ  //liquid water of element 1
//    REAL, INTENT(INOUT) :: WICE  //ice of element 1 [kg/m2]
//    REAL, INTENT(INOUT) :: T     //node temperature of element 1 [k]
//
//// local 
//
//    REAL                :: DZC   //total thickness of nodes 1 and 2 (DZC=DZ+DZ2).
//    REAL                :: WLIQC //combined liquid water [kg/m2]
//    REAL                :: WICEC //combined ice [kg/m2]
			double TC;    //combined node temperature [k]
//    REAL                :: H     //enthalpy of element 1 [J/m2]
//    REAL                :: H2    //enthalpy of element 2 [J/m2]
//    REAL                :: HC    //temporary

//-----------------------------------------------------------------------

			double DZC = DZ + DZ2;
			double WICEC = (WICE + WICE2);
			double WLIQC = (WLIQ + WLIQ2);
			double H = (CICE * WICE + CWAT * WLIQ) * (T - TFRZ) + HFUS * WLIQ;
			double H2 = (CICE * WICE2 + CWAT * WLIQ2) * (T2 - TFRZ) + HFUS * WLIQ2;

			double HC = H + H2;
			if (HC < 0)
				TC = TFRZ + HC / (CICE * WICEC + CWAT * WLIQC);
			else if (HC <= HFUS * WLIQC)
				TC = TFRZ;
			else
				TC = TFRZ + (HC - HFUS * WLIQC) / (CICE * WICEC + CWAT * WLIQC);
    

			DZ = DZC;
			WICE = WICEC;
			WLIQ = WLIQC;
			T = TC;

		}
//  END SUBROUTINE COMBO_GLACIER
		// ==================================================================================================
		public static void DIVIDE_GLACIER(int NSNOW, int NSOIL,                          //in
			ref int ISNOW, FortDoubleArray STC, FortDoubleArray SNICE, FortDoubleArray SNLIQ, FortDoubleArray DZSNSO)  //inout
		{
// ----------------------------------------------------------------------
//    IMPLICIT NONE
// ----------------------------------------------------------------------
// input

//    INTEGER, INTENT(IN)                            :: NSNOW //maximum no. of snow layers [ =3]
//    INTEGER, INTENT(IN)                            :: NSOIL //no. of soil layers [ =4]
//
//// input and output
//
//    INTEGER                        , INTENT(INOUT) :: ISNOW //actual no. of snow layers 
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC   //snow layer temperature [k]
//    REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE //snow layer ice [mm]
//    REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ //snow layer liquid water [mm]
//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO//snow layer depth [m]
//
//// local variables:
//
//    INTEGER                                        :: J     //indices
//    INTEGER                                        :: MSNO  //number of layer (top) to MSNO (bot)
//    REAL                                           :: DRR   //thickness of the combined [m]
//    REAL, DIMENSION(       1:NSNOW)                :: DZ    //snow layer thickness [m]
//    REAL, DIMENSION(       1:NSNOW)                :: SWICE //partial volume of ice [m3/m3]
//    REAL, DIMENSION(       1:NSNOW)                :: SWLIQ //partial volume of liquid water [m3/m3]
//    REAL, DIMENSION(       1:NSNOW)                :: TSNO  //node temperature [k]
			double ZWICE; //temporary
			double ZWLIQ; //temporary
			double PROPOR;//temporary
			double DTDZ;  //temporary
// ----------------------------------------------------------------------

			FortDoubleArray DZ = new FortDoubleArray(1, NSNOW);
			FortDoubleArray SWICE = new FortDoubleArray(1, NSNOW);
			FortDoubleArray SWLIQ = new FortDoubleArray(1, NSNOW);
			FortDoubleArray TSNO = new FortDoubleArray(1, NSNOW);
			//DO J = 1,NSNOW
			for (int J = 1; J <= NSNOW; J++) {
				if (J <= Math.Abs(ISNOW)) {
					DZ[J] = DZSNSO[J + ISNOW];
					SWICE[J] = SNICE[J + ISNOW];
					SWLIQ[J] = SNLIQ[J + ISNOW];
					TSNO[J] = STC[J + ISNOW];
				}
			}

			int MSNO = Math.Abs(ISNOW);

			if (MSNO == 1) {
				// Specify a new snow layer
				if (DZ[1] > 0.05) {
					MSNO = 2;
					DZ[1] = DZ[1] / 2;
					SWICE[1] = SWICE[1] / 2;
					SWLIQ[1] = SWLIQ[1] / 2;
					DZ[2] = DZ[1];
					SWICE[2] = SWICE[1];
					SWLIQ[2] = SWLIQ[1];
					TSNO[2] = TSNO[1];
				}
			}

			if (MSNO > 1) {
				if (DZ[1] > 0.05) {
					double DRR = DZ[1] - 0.05;
					PROPOR = DRR / DZ[1];
					ZWICE = PROPOR * SWICE[1];
					ZWLIQ = PROPOR * SWLIQ[1];
					PROPOR = 0.05 / DZ[1];
					SWICE[1] = PROPOR * SWICE[1];
					SWLIQ[1] = PROPOR * SWLIQ[1];
					DZ[1] = 0.05;
					double dz2 = DZ[2];
					double swliq = SWLIQ[2];
					double swice = SWICE[2];
					double tsno = TSNO[2];
					COMBO_GLACIER(ref dz2, ref swliq, ref swice, ref tsno, DRR, 
						ZWLIQ, ZWICE, TSNO[1]);
					DZ[2] = dz2;
					SWLIQ[2] = swliq;
					SWICE[2] = swice;
					TSNO[2] = tsno;
					// subdivide a new layer
//             if (MSNO <= 2 && DZ[2] > 0.20) THEN  // MB: change limit
					if (MSNO <= 2 && DZ[2] > 0.10) {
						MSNO = 3;
						DTDZ = (TSNO[1] - TSNO[2]) / ((DZ[1] + DZ[2]) / 2);
						DZ[2] = DZ[2] / 2;
						SWICE[2] = SWICE[2] / 2;
						SWLIQ[2] = SWLIQ[2] / 2;
						DZ[3] = DZ[2];
						SWICE[3] = SWICE[2];
						SWLIQ[3] = SWLIQ[2];
						TSNO[3] = TSNO[2] - DTDZ * DZ[2] / 2;
						if (TSNO[3] >= TFRZ)
							TSNO[3] = TSNO[2];
						else
							TSNO[2] = TSNO[2] + DTDZ * DZ[2] / 2;
                

					}
				}
			}

			if (MSNO > 2) {
				if (DZ[2] > 0.2) {
					double DRR = DZ[2] - 0.2;
					PROPOR = DRR / DZ[2];
					ZWICE = PROPOR * SWICE[2];
					ZWLIQ = PROPOR * SWLIQ[2];
					PROPOR = 0.2 / DZ[2];
					SWICE[2] = PROPOR * SWICE[2];
					SWLIQ[2] = PROPOR * SWLIQ[2];
					DZ[2] = 0.2;
					double dz3 = DZ[3];
					double swliq = SWLIQ[3];
					double swice = SWICE[3];
					double tsno = TSNO[3];
					COMBO_GLACIER(ref dz3, ref swliq, ref swice, ref tsno, DRR, ZWLIQ, ZWICE, TSNO[2]);
					DZ[3] = dz3;
					SWLIQ[3] = swliq;
					SWICE[3] = swice;
					TSNO[3] = tsno;
				}
			}

			ISNOW = -MSNO;

			//DO J = ISNOW+1,0
			for (int J = ISNOW + 1; J <= 0; J++) {
				DZSNSO[J] = DZ[J - ISNOW];
				SNICE[J] = SWICE[J - ISNOW];
				SNLIQ[J] = SWLIQ[J - ISNOW];
				STC[J] = TSNO[J - ISNOW];
			}


//    DO J = ISNOW+1,NSOIL
//    WRITE(*,"(I5,7F10.3)") J, DZSNSO[J], SNICE[J], SNLIQ[J],STC[J]
//    END DO

		}
//  END SUBROUTINE DIVIDE_GLACIER
		// ==================================================================================================
		public static void SNOWH2O_GLACIER(int NSNOW, int NSOIL, double DT, double QSNFRO, double QSNSUB,  //in 
			double QRAIN,                                  //in
			ref int ISNOW, FortDoubleArray DZSNSO, ref double SNOWH, ref double SNEQV, FortDoubleArray SNICE,  //inout
			FortDoubleArray SNLIQ, FortDoubleArray SH2O, FortDoubleArray SICE, FortDoubleArray STC,          //inout
			ref double PONDING1, ref double PONDING2,          //inout
			out double QSNBOT)                                   //out
		{
// ----------------------------------------------------------------------
// Renew the mass of ice lens (SNICE) and liquid (SNLIQ) of the
// surface snow layer resulting from sublimation (frost) / evaporation (dew)
// ----------------------------------------------------------------------
   
// ----------------------------------------------------------------------
// input

//   INTEGER,                         INTENT(IN)    :: NSNOW  //maximum no. of snow layers[=3]
//   INTEGER,                         INTENT(IN)    :: NSOIL  //No. of soil layers[=4]
//   REAL,                            INTENT(IN)    :: DT     //time step
//   REAL,                            INTENT(IN)    :: QSNFRO //snow surface frost rate[mm/s]
//   REAL,                            INTENT(IN)    :: QSNSUB //snow surface sublimation rate[mm/s]
//   REAL,                            INTENT(IN)    :: QRAIN  //snow surface rain rate[mm/s]
//
//// output
//
//   REAL,                            INTENT(OUT)   :: QSNBOT //melting water out of snow bottom [mm/s]
//
//// input and output
//
//   INTEGER,                         INTENT(INOUT) :: ISNOW  //actual no. of snow layers
//   REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO // snow layer depth [m]
//   REAL,                            INTENT(INOUT) :: SNOWH  //snow height [m]
//   REAL,                            INTENT(INOUT) :: SNEQV  //snow water eqv. [mm]
//   REAL, DIMENSION(-NSNOW+1:0),     INTENT(INOUT) :: SNICE  //snow layer ice [mm]
//   REAL, DIMENSION(-NSNOW+1:0),     INTENT(INOUT) :: SNLIQ  //snow layer liquid water [mm]
//   REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O   //soil liquid moisture (m3/m3)
//   REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SICE   //soil ice moisture (m3/m3)
//   REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC    //snow layer temperature [k]
//   REAL,                            INTENT(INOUT) :: PONDING1
//   REAL,                            INTENT(INOUT) :: PONDING2
//
//// local variables:
//
//   INTEGER                     :: J         //do loop/array indices
//   REAL                        :: QIN       //water flow into the element (mm/s)
//   REAL                        :: QOUT      //water flow out of the element (mm/s)
//   REAL                        :: WGDIF     //ice mass after minus sublimation
//   REAL, DIMENSION(-NSNOW+1:0) :: VOL_LIQ   //partial volume of liquid water in layer
//   REAL, DIMENSION(-NSNOW+1:0) :: VOL_ICE   //partial volume of ice lens in layer
//   REAL, DIMENSION(-NSNOW+1:0) :: EPORE     //effective porosity = porosity - VOL_ICE
//   REAL :: PROPOR, TEMP
// ----------------------------------------------------------------------
			FortDoubleArray VOL_LIQ = new FortDoubleArray(-NSNOW + 1, 0);
			FortDoubleArray VOL_ICE = new FortDoubleArray(-NSNOW + 1, 0);
			FortDoubleArray EPORE = new FortDoubleArray(-NSNOW + 1, 0);
//for the case when SNEQV becomes "0" after "COMBINE"

			if (SNEQV == 0) {
				SICE[1] = SICE[1] + (QSNFRO - QSNSUB) * DT / (DZSNSO[1] * 1000);
			}

// for shallow snow without a layer
// snow surface sublimation may be larger than existing snow mass. To conserve water,
// excessive sublimation is used to reduce soil water. Smaller time steps would tend 
// to aviod this problem.

			if (ISNOW == 0 && SNEQV > 0) {
				double TEMP = SNEQV;
				SNEQV = SNEQV - QSNSUB * DT + QSNFRO * DT;
				double PROPOR = SNEQV / TEMP;
				SNOWH = Math.Max(0, PROPOR * SNOWH);

				if (SNEQV < 0) {
					SICE[1] = SICE[1] + SNEQV / (DZSNSO[1] * 1000);
					SNEQV = 0;
					SNOWH = 0;
				}
				if (SICE[1] < 0) {
					SH2O[1] = SH2O[1] + SICE[1];
					SICE[1] = 0;
				}
			}

			if (SNOWH <= 1E-8 || SNEQV <= 1E-6) {
				SNOWH = 0.0;
				SNEQV = 0.0;
			}

// for deep snow

			if (ISNOW < 0) { //KWM added this if statement to prevent out-of-bounds array references

				double WGDIF = SNICE[ISNOW + 1] - QSNSUB * DT + QSNFRO * DT;
				SNICE[ISNOW + 1] = WGDIF;
				if (WGDIF < 1e-6 && ISNOW < 0) {
					COMBINE_GLACIER(NSNOW, NSOIL,                          //in
						ref ISNOW, SH2O, STC, SNICE, SNLIQ,  //inout
						DZSNSO, SICE, ref SNOWH, ref SNEQV,          //inout
						ref PONDING1, ref PONDING2);                        //inout
				}
				//KWM:  Subroutine COMBINE can change ISNOW to make it 0 again?
				if (ISNOW < 0) { //KWM added this if statement to prevent out-of-bounds array references
					SNLIQ[ISNOW + 1] = SNLIQ[ISNOW + 1] + QRAIN * DT;
					SNLIQ[ISNOW + 1] = Math.Max(0, SNLIQ[ISNOW + 1]);
				}
      
			} //KWM  -- Can the ENDIF be moved toward the end of the subroutine (Just set QSNBOT=0)?

// Porosity and partial volume

			//KWM Looks to me like loop index / if test can be simplified.

   
			//DO J = -NSNOW+1, 0
			for (int J = -NSNOW + 1; J <= 0; J++) {
				if (J >= ISNOW + 1) {
					VOL_ICE[J] = Math.Min(1, SNICE[J] / (DZSNSO[J] * DENICE));
					EPORE[J] = 1 - VOL_ICE[J];
					VOL_LIQ[J] = Math.Min(EPORE[J], SNLIQ[J] / (DZSNSO[J] * DENH2O));
				}
			}

			double QIN = 0;
			double QOUT = 0;

			//KWM Looks to me like loop index / if test can be simplified.

			for (int J = -NSNOW + 1; J <= 0; J++) {
   	
				if (J >= ISNOW + 1) {
					SNLIQ[J] += QIN;
					if (J <= -1) {
						if (EPORE[J] < 0.05 || EPORE[J + 1] < 0.05) {
							QOUT = 0;
						} else {
							QOUT = Math.Max(0, (VOL_LIQ[J] - SSI * EPORE[J]) * DZSNSO[J]);
							QOUT = Math.Min(QOUT, (1 - VOL_ICE[J + 1] - VOL_LIQ[J + 1]) * DZSNSO[J + 1]);
						}
					} else {
						QOUT = Math.Max(0, (VOL_LIQ[J] - SSI * EPORE[J]) * DZSNSO[J]);
					}
					QOUT = QOUT * 1000;
					SNLIQ[J] = SNLIQ[J] - QOUT;
					QIN = QOUT;
				}
			}

// Liquid water from snow bottom to soil

			QSNBOT = QOUT / DT;           // mm/s

		}
//  END SUBROUTINE SNOWH2O_GLACIER
		// ********************* end of water subroutines ******************************************
		// ==================================================================================================
		public static void ERROR_GLACIER(int ILOC, int JLOC, double SWDOWN, double FSA, double FSR, double FIRA, 
			double FSH, double FGEV, double SSOIL, double SAG, double PRCP, double EDIR, 
			double  RUNSRF, double RUNSUB, double SNEQV, double DT, double BEG_WB)
		{
// --------------------------------------------------------------------------------------------------
// check surface energy balance and water balance
// --------------------------------------------------------------------------------------------------
//  IMPLICIT NONE
// --------------------------------------------------------------------------------------------------
// inputs
//  INTEGER                        , INTENT(IN) :: ILOC   //grid index
//  INTEGER                        , INTENT(IN) :: JLOC   //grid index
//  REAL                           , INTENT(IN) :: SWDOWN //downward solar filtered by sun angle [w/m2]
//  REAL                           , INTENT(IN) :: FSA    //total absorbed solar radiation (w/m2)
//  REAL                           , INTENT(IN) :: FSR    //total reflected solar radiation (w/m2)
//  REAL                           , INTENT(IN) :: FIRA   //total net longwave rad (w/m2)  [+ to atm]
//  REAL                           , INTENT(IN) :: FSH    //total sensible heat (w/m2)     [+ to atm]
//  REAL                           , INTENT(IN) :: FGEV   //ground evaporation heat (w/m2) [+ to atm]
//  REAL                           , INTENT(IN) :: SSOIL  //ground heat flux (w/m2)        [+ to soil]
//  REAL                           , INTENT(IN) :: SAG
//
//  REAL                           , INTENT(IN) :: PRCP   //precipitation rate (kg m-2 s-1)
//  REAL                           , INTENT(IN) :: EDIR   //soil surface evaporation rate[mm/s]
//  REAL                           , INTENT(IN) :: RUNSRF //surface runoff [mm/s] 
//  REAL                           , INTENT(IN) :: RUNSUB //baseflow (saturation excess) [mm/s]
//  REAL                           , INTENT(IN) :: SNEQV  //snow water eqv. [mm]
//  REAL                           , INTENT(IN) :: DT     //time step [sec]
//  REAL                           , INTENT(IN) :: BEG_WB //water storage at begin of a timesetp [mm]
//
//  REAL                                        :: END_WB //water storage at end of a timestep [mm]
//  REAL                                        :: ERRWAT //error in water balance [mm/timestep]
//  REAL                                        :: ERRENG //error in surface energy balance [w/m2]
//  REAL                                        :: ERRSW  //error in shortwave radiation balance [w/m2]
//  CHARACTER(len=256)                          :: message
// --------------------------------------------------------------------------------------------------
			double ERRSW = SWDOWN - (FSA + FSR);
			if (ERRSW > 0.01) {          // w/m2
//     WRITE(*,*) "SAG    =",SAG
//     WRITE(*,*) "FSA    =",FSA
//     WRITE(*,*) "FSR    =",FSR
//     WRITE(message,*) "ERRSW =",ERRSW
//     call wrf_message(trim(message))
//     call wrf_error_fatal("Radiation budget problem in NOAHMP GLACIER")
			}

			double ERRENG = SAG - (FIRA + FSH + FGEV + SSOIL);
			if (ERRENG > 0.01) {
//      write(message,*) "ERRENG =",ERRENG
//      call wrf_message(trim(message))
//      WRITE(message,"(i6,1x,i6,1x,5F10.4)")ILOC,JLOC,SAG,FIRA,FSH,FGEV,SSOIL
//      call wrf_message(trim(message))
//      call wrf_error_fatal("Energy budget problem in NOAHMP GLACIER")
			}

			double END_WB = SNEQV;
			double ERRWAT = END_WB - BEG_WB - (PRCP - EDIR - RUNSRF - RUNSUB) * DT;

//ifndef WRF_HYDRO
			if (Math.Abs(ERRWAT) > 0.1) {
				if (ERRWAT > 0) {
					Console.WriteLine("The model is gaining water (ERRWAT is positive)");
				} else {
					Console.WriteLine("The model is losing water (ERRWAT is negative)");
				}
//      write(message, *) "ERRWAT =",ERRWAT, "kg m{-2} timestep{-1}"
//      call wrf_message(trim(message))
//      WRITE(message,"("    I      J     END_WB     BEG_WB       PRCP       EDIR      RUNSRF     RUNSUB")")
//           call wrf_message(trim(message))
//           WRITE(message,"(i6,1x,i6,1x,2f15.3,4f11.5)")ILOC,JLOC,END_WB,BEG_WB,PRCP*DT,
//                EDIR*DT,RUNSRF*DT,RUNSUB*DT
//           call wrf_message(trim(message))
//           call wrf_error_fatal("Water budget problem in NOAHMP GLACIER")
//        END if
//endif

			}// END SUBROUTINE ERROR_GLACIER
// ==================================================================================================

//  public static void NOAHMP_OPTIONS_GLACIER(idveg     ,iopt_crs  ,iopt_btr  ,iopt_run  ,iopt_sfc  ,iopt_frz ,  
//                             iopt_inf  ,iopt_rad  ,iopt_alb  ,iopt_snf  ,iopt_tbot, iopt_stc )
//
////  IMPLICIT NONE
////
////  INTEGER,  INTENT(IN) :: idveg     //dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1
////  INTEGER,  INTENT(IN) :: iopt_crs  //canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
////  INTEGER,  INTENT(IN) :: iopt_btr  //soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
////  INTEGER,  INTENT(IN) :: iopt_run  //runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
////  INTEGER,  INTENT(IN) :: iopt_sfc  //surface layer drag coeff (CH  CM) (1->M-O; 2->Chen97)
////  INTEGER,  INTENT(IN) :: iopt_frz  //supercooled liquid water (1-> NY06; 2->Koren99)
////  INTEGER,  INTENT(IN) :: iopt_inf  //frozen soil permeability (1-> NY06; 2->Koren99)
////  INTEGER,  INTENT(IN) :: iopt_rad  //radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
////  INTEGER,  INTENT(IN) :: iopt_alb  //snow surface albedo (1->BATS; 2->CLASS)
////  INTEGER,  INTENT(IN) :: iopt_snf  //rainfall  snowfall (1-Jordan91; 2->BATS; 3->Noah)
////  INTEGER,  INTENT(IN) :: iopt_tbot //lower boundary of soil temperature (1->zero-flux; 2->Noah)
////
////  INTEGER,  INTENT(IN) :: iopt_stc  //snow/soil temperature time scheme (only layer 1)
//                                    // 1 -> semi-implicit; 2 -> full implicit (original Noah)
//
//// -------------------------------------------------------------------------------------------------
//
//  dveg = idveg
//  
//  opt_crs  = iopt_crs  
//  opt_btr  = iopt_btr  
//  opt_run  = iopt_run  
//  opt_sfc  = iopt_sfc  
//  opt_frz  = iopt_frz  
//  opt_inf  = iopt_inf  
//  opt_rad  = iopt_rad  
//  opt_alb  = iopt_alb  
//  opt_snf  = iopt_snf  
//  opt_tbot = iopt_tbot 
//  opt_stc  = iopt_stc
//  
//  end subroutine noahmp_options_glacier
// 
//
//// ==================================================================================================
//
//
//	}
		}
	}
}


