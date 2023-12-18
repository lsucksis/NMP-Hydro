/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2018/1/10
 * Time: 10:26
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;

namespace NoahMP
{
	public class ENERGY_Input
	{
		// inputs
		public int ILOC;
		public int JLOC;
		public int ICE;
		//ice (ice = 1)
		
		/// <summary>
		/// vegetation physiology type
		/// </summary>
		public int VEGTYP;
		/// <summary>
		/// surface type: 1->soil; 2->lake
		/// </summary>
		public int IST;
		/// <summary>
		/// soil color type (1-lighest; 8-darkest)
		/// </summary>
		public int ISC;
		/// <summary>
		/// maximum no. of snow layers
		/// </summary>
		public static int NSNOW;
		/// <summary>
		/// number of soil layers
		/// </summary>
		public static int NSOIL;
		/// <summary>
		/// number of root layers
		/// </summary>
		public static int NROOT;
		/// <summary>
		/// actual no. of snow layers
		/// </summary>
		public int ISNOW;
		/// <summary>
		/// time step [sec]
		/// </summary>
		public double DT;
		/// <summary>
		/// snowfall on the ground (mm/s)
		/// </summary>
		public double QSNOW;
		/// <summary>
		/// density air (kg/m3)
		/// </summary>
		public double RHOAIR;
		/// <summary>
		/// vapor pressure air (pa)
		/// </summary>
		public double EAIR;
		/// <summary>
		/// pressure (pa)
		/// </summary>
		public double SFCPRS;
		/// <summary>
		/// specific humidity (kg/kg)
		/// </summary>
		public double QAIR;
		/// <summary>
		/// air temperature [K]
		/// </summary>
		public double SFCTMP;
		/// <summary>
		/// potential temperature [K]
		/// </summary>
		public double THAIR;
		/// <summary>
		/// downward longwave radiation (w/m2)
		/// </summary>
		public double LWDN;
		/// <summary>
		/// wind speed in e-w dir (m/s)
		/// </summary>
		public double UU;
		/// <summary>
		/// wind speed in n-s dir (m/s)
		/// </summary>
		public double VV;
		/// <summary>
		/// incoming direct solar rad. (w/m2)
		/// </summary>
		public double[] SOLAD = new double[2];
		/// <summary>
		/// incoming diffuse solar rad. (w/m2)
		/// </summary>
		public double[] SOLAI = new double[2];
		/// <summary>
		/// cosine solar zenith angle (0-1)
		/// </summary>
		public double COSZ;
		/// <summary>
		/// LAI adjusted for burying by snow
		/// </summary>
		public double ELAI;
		/// <summary>
		/// LAI adjusted for burying by snow
		/// </summary>
		public double ESAI;
		/// <summary>
		/// vol. soil heat capacity [j/m3/k]
		/// </summary>
		public double CSOIL;
		/// <summary>
		/// fraction of canopy that is wet [-]
		/// </summary>
		public double FWET;
		/// <summary>
		/// top of canopy layer (m)
		/// </summary>
		public double HTOP;
		/// <summary>
		/// greeness vegetation fraction (-)
		/// </summary>
		public double FVEG;
		/// <summary>
		/// latitude (radians)
		/// </summary>
		public double LAT;
		/// <summary>
		/// canopy-intercepted liquid water (mm)
		/// </summary>
		public double CANLIQ;
		/// <summary>
		/// canopy-intercepted ice mass (mm)
		/// </summary>
		public double CANICE;
		/// <summary>
		/// foliage nitrogen (%)
		/// </summary>
		public double FOLN;
		/// <summary>
		/// atmospheric co2 concentration (pa)
		/// </summary>
		public double CO2AIR;
		/// <summary>
		/// atmospheric o2 concentration (pa)
		/// </summary>
		public double O2AIR;
		/// <summary>
		/// growing season index (0=off, 1=on)
		/// </summary>
		public double IGS;
		/// <summary>
		/// reference height (m)
		/// </summary>

		public double ZREF;
		/// <summary>
		/// bottom condition for soil temp. [K]
		/// </summary>
		public double TBOT;
		/// <summary>
		/// depth for TBOT [m]
		/// </summary>
		public double ZBOT;
		/// <summary>
		/// layer-bottom depth from snow surf [m], DIMENSION(-NSNOW+1:NSOIL),
		/// </summary>
		public double[] ZSNSO = new double[NSNOW + NSOIL];
		/// <summary>
		/// layer-bottom depth from soil surf [m], DIMENSION(       1:NSOIL),
		/// </summary>
		public double[] ZSOIL = new double[NSOIL];
		/// <summary>
		/// depth of snow & soil layer-bottom [m], DIMENSION(-NSNOW+1:NSOIL),
		/// </summary>
		public double[] DZSNSO = new double[NSNOW + NSOIL];
		

		//jref:start; in
		public int ISURBAN;
		public int IZ0TLND;
		/// <summary>
		/// cloud water mixing ratio
		/// </summary>
		public double QC;
		/// <summary>
		/// planetary boundary layer height
		/// </summary>
		public double PBLH;
		/// <summary>
		/// mixing ratio at lowest model layer
		/// </summary>
		public double QSFC;
		/// <summary>
		/// pressure at lowest model layer
		/// </summary>
		public double PSFC;
		/// <summary>
		/// horisontal resolution
		/// </summary>
		public double DX;
		/// <summary>
		/// thickness of lowest layer
		/// </summary>
		public double DZ8W;
		/// <summary>
		/// mixing ratio (kg/kg)
		/// </summary>
		public double Q2;
		

		// outputs
		/// <summary>
		/// phase change index [1-melt; 2-freeze], DIMENSION(-NSNOW+1:NSOIL)
		/// </summary>
		public int[] IMELT = new int[NSNOW + NSOIL];
		/// <summary>
		/// partial volume ice [m3/m3], DIMENSION(-NSNOW+1:    0)
		/// </summary>
		public double[] SNICEV = new double[NSNOW];
		/// <summary>
		/// partial volume liq. water [m3/m3] , DIMENSION(-NSNOW+1:    0)
		/// </summary>
		public double[] SNLIQV = new double[NSNOW];
		/// <summary>
		/// effective porosity [m3/m3], DIMENSION(-NSNOW+1:    0)
		/// </summary>
		public double[] EPORE = new double[NSNOW];
		/// <summary>
		/// snow cover fraction (-)
		/// </summary>
		public double FSNO;
		/// <summary>
		/// snowmelt [mm/s]
		/// </summary>
		public double QMELT;
		/// <summary>
		/// pounding at ground [mm]
		/// </summary>
		public double PONDING;
		/// <summary>
		/// solar rad. absorbed by veg. (w/m2)
		/// </summary>
		public double SAV;
		/// <summary>
		/// solar rad. absorbed by ground (w/m2)
		/// </summary>
		public double SAG;
		/// <summary>
		/// tot. absorbed solar radiation (w/m2)
		/// </summary>
		public double FSA;
		/// <summary>
		/// tot. reflected solar radiation (w/m2)
		/// </summary>
		public double FSR;
		/// <summary>
		/// wind stress: e-w (n/m2)
		/// </summary>
		public double TAUX;
		/// <summary>
		/// wind stress: n-s (n/m2)
		/// </summary>
		public double TAUY;
		/// <summary>
		/// total net LW. rad (w/m2)   [+ to atm]
		/// </summary>
		public double FIRA;
		/// <summary>
		/// total sensible heat (w/m2) [+ to atm]
		/// </summary>
		public double FSH;
		/// <summary>
		/// canopy evaporation (w/m2)  [+ to atm]
		/// </summary>
		public double FCEV;
		/// <summary>
		/// ground evaporation (w/m2)  [+ to atm]
		/// </summary>
		public double FGEV;
		/// <summary>
		/// transpiration (w/m2)       [+ to atm]
		/// </summary>
		public double FCTR;
		/// <summary>
		/// radiative temperature [K]
		/// </summary>
		public double TRAD;
		/// <summary>
		/// 2 m height air temperature [K]
		/// </summary>
		public double T2M;
		/// <summary>
		/// total photosyn. (umolco2/m2/s) [+]
		/// </summary>
		public double PSN;
		/// <summary>
		/// total photosyn. active energy (w/m2)
		/// </summary>
		public double APAR;
		/// <summary>
		/// ground heat flux (w/m2)   [+ to soil]
		/// </summary>
		public double SSOIL;
		/// <summary>
		/// soil water transpiration factor (0-1), DIMENSION(       1:NSOIL)
		/// </summary>
		public double[] BTRANI = new double[NSNOW];
		/// <summary>
		/// soil water transpiration factor (0-1)
		/// </summary>
		public double BTRAN;
		//
		//  public double                                  LATHEA //latent heat vap./sublimation (j/kg)
		/// <summary>
		/// latent heat vap./sublimation (j/kg)
		/// </summary>
		public double LATHEAV;
		/// <summary>
		/// latent heat vap./sublimation (j/kg)
		/// </summary>
		public double LATHEAG;
		/// <summary>
		/// used to define latent heat pathway
		/// </summary>
		bool FROZEN_GROUND;
		/// <summary>
		/// used to define latent heat pathway
		/// </summary>
		bool FROZEN_CANOPY;
		
		/// <summary>
		/// veg. reflected solar radiation (w/m2)
		/// </summary>
		public double FSRV;
		/// <summary>
		/// ground reflected solar radiation (w/m2)
		/// </summary>
		public double FSRG;
		/// <summary>
		/// sunlit leaf stomatal resistance (s/m)
		/// </summary>
		public double RSSUN;
		/// <summary>
		/// shaded leaf stomatal resistance (s/m)
		/// </summary>
		public double RSSHA;
		/// <summary>
		/// 2-m air temperature over vegetated part [k]
		/// </summary>
		public double T2MV;
		/// <summary>
		/// 2-m air temperature over bare ground part [k]
		/// </summary>
		public double T2MB;
		/// <summary>
		/// 
		/// </summary>
		public double BGAP;
		public double WGAP;
		
		
		public double TS;
		//surface temperature [K]
		public double TV;
		//vegetation temperature [K]
		public double TG;
		//ground temperature [K]
		public double[] STC = new double[NSNOW + NSOIL];
		//snow/soil temperature [k], DIMENSION(-NSNOW+1:NSOIL)
		public double SNOWH;
		//snow height [m]
		public double SNEQV;
		//snow mass (mm)
		public double SNEQVO;
		//snow mass at last time step (mm)
		public double[] SH2O = new double[NSOIL];
		//liquid soil moisture [m3/m3] , DIMENSION(       1:NSOIL)
		public double[] SMC = new double[NSOIL];
		//soil moisture (ice + liq.) [m3/m3], DIMENSION(       1:NSOIL)
		public double[] SNICE = new double[NSNOW];
		//snow ice mass (kg/m2), DIMENSION(-NSNOW+1:    0)
		public double[] SNLIQ = new double[NSNOW];
		//snow liq mass (kg/m2), DIMENSION(-NSNOW+1:    0)
		public double EAH;
		//canopy air vapor pressure (pa)
		public double TAH;
		//canopy air temperature [K]
		public double ALBOLD;
		//snow albedo at last time step(CLASS type)
		public double TAUSS;
		//non-dimensional snow age
		public double CM;
		//momentum drag coefficient
		public double CH;
		//sensible heat exchange coefficient
		public double Q1;

		public double EMISSI;

	}
}
