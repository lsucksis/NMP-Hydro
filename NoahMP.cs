/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2016/11/26
 * Time: 20:58
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 * Completed in Jan 2018
 */
using System;

namespace NoahMP
{
	/// <summary>
	/// NoahMP的原始类，用于初始化里面的变量。其继承类里面主要是方法
	/// </summary>
	public class NoahMP
	{
		/// <summary>
		/// maximum no. of snow layers
		/// </summary>
		public static int NSNOW;
		/// <summary>
		/// number of soil layers
		/// </summary>
		public static int NSOIL;


		VEG_PARAMS vegparams = new VEG_PARAMS();
		/// <summary>
		/// acceleration due to gravity (m/s2)
		/// </summary>
		public static  double GRAV = 9.80616;
		/// <summary>
		/// Stefan-Boltzmann constant (w/m2/k4)
		/// </summary>
		public static  double SB = 5.67E-08;
		/// <summary>
		/// von Karman constant
		/// </summary>
		public static  double VKC = 0.40;
		/// <summary>
		/// freezing/melting point [K]
		/// </summary>
		public static  double TFRZ = 273.16;
		/// <summary>
		/// latent heat of sublimation (j/kg)
		/// </summary>
		public static  double HSUB = 2.8440E06;
		/// <summary>
		/// latent heat of vaporization (j/kg)
		/// </summary>
		public static  double HVAP = 2.5104E06;
		/// <summary>
		/// latent heat of fusion (j/kg)
		/// </summary>
		public static double HFUS = 0.3336E06;
		/// <summary>
		/// specific heat capacity of water (j/m3/k)
		/// </summary>
		public static double CWAT = 4.188E06;
		/// <summary>
		/// specific heat capacity of ice (j/m3/k)
		/// </summary>
		public static  double CICE = 2.094E06;
		/// <summary>
		/// heat capacity dry air at const pres (j/kg/k)
		/// </summary>
		public static double CPAIR = 1004.64;
		/// <summary>
		/// thermal conductivity of water (w/m/k)
		/// </summary>
		public static  double TKWAT = 0.6;
		/// <summary>
		/// thermal conductivity of ice (w/m/k)
		/// </summary>
		public static double TKICE = 2.2;
		/// <summary>
		/// thermal conductivity of air (w/m/k)
		/// </summary>
		public static double TKAIR = 0.023;
		/// <summary>
		/// gas constant for dry air (j/kg/k)
		/// </summary>
		public static double RAIR = 287.04;
		/// <summary>
		/// gas constant for  water vapor (j/kg/k)
		/// </summary>
		public static double RW = 461.269;
		/// <summary>
		/// density of water (kg/m3)
		/// </summary>
		public static double DENH2O = 1000;
		/// <summary>
		/// density of ice (kg/m3)
		/// </summary>
		public static double DENICE = 917;
		//
		
		// runoff parameters used for SIMTOP and SIMGM:
		/// <summary>
		/// gridcell mean topgraphic index (global mean)
		/// </summary>
		public static	double TIMEAN = 10.5;
		/// <summary>
		/// maximum surface saturated fraction (global mean)
		/// </summary>
		public static	double FSATMX = 0.38;
		

		// adjustable parameters for snow processes
		/// <summary>
		/// melting factor (-)
		/// </summary>
		public static double M = 2.50;
		/// <summary>
		/// snow surface roughness length (m) (0.002)
		/// </summary>
		public static double Z0SNO = 0.002;
		/// <summary>
		/// liquid water holding capacity for snowpack (m3/m3) (0.03)
		/// </summary>
		public static  double SSI = 0.03;
		/// <summary>
		/// new snow mass to fully cover old snow (mm) equivalent to 10mm depth (density = 100 kg/m3)
		/// </summary>
		public static  double SWEMX = 1.00;
		//
		//

		//public  int NROOT = 0;
		
		/// <summary>
		/// 土壤类型数量
		/// </summary>		
		public static int NSLTYPE = 30;
		/// <summary>
		/// vertical crown radius (m),B parameter ( Assigned in REDPRM )
		/// </summary>
		public static double[]	BB = new double[NSLTYPE];
		/// <summary>
		/// dry soil moisture threshold where direct evap from top
		/// layer ends (volumetric) ( Assigned in REDPRM )
		/// </summary>
		public static  double[]	DRYSMC = new double[NSLTYPE];
		/// <summary>
		/// soil thermal diffusivity/conductivity coef ( Assigned in REDPRM )
		/// </summary>
		public static double[]	F11 = new double[NSLTYPE];
		/// <summary>
		/// 应该是饱和土壤含水率(体积比)
		/// porosity, saturated value of soil moisture (volumetric)
		/// </summary>
		public static double[]	MAXSMC = new double[NSLTYPE];
		/// <summary>
		/// reference soil moisture (field capacity) (volumetric) ( Assigned in REDPRM )
		/// 参考土壤含水率(体积比)
		/// </summary>
		public static double[]	REFSMC = new double[NSLTYPE];
		/// <summary>
		/// saturated soil matric potential ( Assigned in REDPRM )
		/// </summary>
		public static double[]	SATPSI = new double[NSLTYPE];
		/// <summary>
		/// saturated soil hydraulic conductivity ( Assigned in REDPRM )
		/// </summary>
		public static  double[]	SATDK = new double[NSLTYPE];
		/// <summary>
		/// saturated soil hydraulic diffusivity ( Assigned in REDPRM )
		/// </summary>
		public static double[]	SATDW = new double[NSLTYPE];
		/// <summary>
		/// wilting point soil moisture (volumetric) ( Assigned in REDPRM )
		/// </summary>
		public static double[]	WLTSMC = new double[NSLTYPE];
		/// <summary>
		/// soil quartz content ( Assigned in REDPRM )
		/// </summary>
		public static double[]	QTZ = new double[NSLTYPE];
		
		//	 REAL, DIMENSION (1:NSLTYPE) :: BB,DRYSMC,F11,                           &
		//        MAXSMC, REFSMC,SATPSI,SATDK,SATDW, WLTSMC,QTZ
		//------------------------------------------------------------------------------------------//
		// From the GENPARM.TBL file
		//------------------------------------------------------------------------------------------//
		/// <summary>
		/// slope index (0 - 1) ( Assigned in REDPRM )
		/// </summary>
		//public 	double SLOPE = 0.1;
		/// <summary>
		/// vol. soil heat capacity [j/m3/K] ( Assigned in REDPRM )
		/// </summary>
		public static  double CSOIL = 2e6;
		/// <summary>
		/// Depth (m) of lower boundary soil temperature ( Assigned in REDPRM )
		/// </summary>
		public static  double ZBOT = -8;
		/// <summary>
		/// Calculate roughness length of heat ( Assigned in REDPRM )
		/// </summary>
		public static  double CZIL = 0.1;
		/// <summary>
		/// used in compute maximum infiltration rate (in INFIL) ( Assigned in REDPRM )
		/// </summary>
		public  static double REFDK = 2e-6;
		/// <summary>
		/// 被Arnault et al. (2015) and Yucel et al. (2015)认为是影响地表径流最重要的参数
		/// </summary>
		public static 	double REFKDT = 1.0;
		/// <summary>
		/// used in compute maximum infiltration rate (in INFIL) ( Assigned in REDPRM )
		/// </summary>
		public static	double FRZK = 0.15;
		//


		// =====================================options for different schemes================================
		

		/// <summary>
		/// 动态植被选项
		/// 1 -> off (use table LAI; use FVEG = SHDFAC from input)
		/// 2 -> on (together with OPT_CRS = 1)
		/// 3 -> off (use table LAI; calculate FVEG)
		/// 4 -> off (use table LAI; use maximum vegetation fraction)
		/// 原NoahMP建议选4
		/// </summary>
		public static 	int DVEG = 1;
	

		/// <summary>
		/// options for canopy stomatal resistance
		///  1-> Ball-Berry; 2->Jarvis
		/// 默认建议选1.选项1和2的计算都已测试过。曾在选项2中修正了两个错误
		/// </summary>
		public static 	int OPT_CRS = 1;
		//= 1    //(must 1 when DVEG = 2)



		/// <summary>
		/// options for soil moisture factor for stomatal resistance
		/// 1-> Noah (soil moisture)
		/// 2-> CLM  (matric potential)
		/// 3-> SSiB (matric potential)
		/// (suggested 1)
		/// 选项1应该已默认测试过了，2和3需要测试一下，但比较简单，出错机率小
		/// </summary>
		public static int OPT_BTR = 1;
		
		
		
		/// <summary>
		///  options for runoff and groundwater
		/// 1 -> TOPMODEL with groundwater, SIMGM (Niu et al. 2007 JGR) ;
		/// 2 -> TOPMODEL with an equilibrium water table (Niu et al. 2005 JGR) ;
		/// 3 -> original surface and subsurface runoff (free drainage)
		/// 4 -> BATS surface and subsurface runoff (free drainage)
		/// 5 -> Miguez-Macho&Fan groundwater scheme (Miguez-Macho et al. 2007 JGR, lateral flow: Fan et al. 2007 JGR)
		/// 其中1，2，3，4已测试通过，5太复杂，未测试
		/// (suggested 1)
		/// </summary>
		public static  int OPT_RUN = 4;
		
		/// <summary>
		/// options for surface layer drag coeff (CH & CM)
		/// 1->M-O ; 2->original Noah (Chen97); 3->MYJ consistent; 4->YSU consistent.
		/// 选项1，2和4应该已问题不大，但是3不收敛，不过解决意义已不大
		/// </summary>
		public static  int OPT_SFC = 1;
		

		/// <summary>
		/// options for supercooled liquid water (or ice fraction)
		/// 1-> no iteration (Niu and Yang, 2006 JHM); 2: Koren's iteration
		/// 1和2均已验证过，问题不大
		/// </summary>
		public static int OPT_FRZ = 2;
		//= 1    //(1 or 2)

		/// <summary>
		/// 冻土渗透选项
		/// options for frozen soil permeability
		/// 1 -> linear effects, more permeable (Niu and Yang, 2006, JHM)
		/// 2 -> nonlinear effects, less permeable (old)
		/// 最先测试的是1，即默认选项。选项2也确认没有问题了
		/// </summary>
		public static  int OPT_INF = 2;
		//= 1    //(suggested 1)

		/// <summary>
		/// options for radiation transfer
		/// 1 -> modified two-stream (gap = F(solar angle, 3D structure ...)<1-FVEG)
		/// 2 -> two-stream applied to grid-cell (gap = 0)
		/// 3 -> two-stream applied to vegetated fraction (gap=1-FVEG)
		/// 选项1已验证过，2和3很简单无需验证
		/// </summary>
		public static int OPT_RAD = 1;
		
		/// <summary>
		/// options for ground snow surface albedo
		/// 1-> BATS; 2 -> CLASS
		/// 两个选项都比较简单，应该不用验证了
		/// </summary>
		public static int OPT_ALB = 2;
		

		/// <summary>
		/// options for partitioning  precipitation into rainfall & snowfall
		/// 1 -> Jordan (1991); 2 -> BATS: when SFCTMP<TFRZ+2.2; 3 -> SFCTMP<TFRZ
		/// 此选项下非常简单，不需要验证
		/// </summary>
		public static  int OPT_SNF = 1;
		//= 1    //(suggested 1)

		/// <summary>
		/// options for lower boundary condition of soil temperature
		/// 1 -> zero heat flux from bottom (ZBOT and TBOT not used)
		/// 2 -> TBOT at ZBOT (8m) read from a file (original Noah)
		/// 该选项下比较简单，已对照过，只要有一个选项下是正确的，则另一个下面也没问题
		/// </summary>
		public static  int OPT_TBOT = 2;
		//= 2   //(suggested 2)

		/// <summary>
		/// options for snow/soil temperature time scheme (only layer 1)
		/// 1 -> semi-implicit; 2 -> full implicit (original Noah)
		/// 此选项下面相对比较简单，应该问题不大
		/// </summary>
		public static  int OPT_STC = 2;
		//= 1    //(suggested 1)
		// ==================================================================================================

		//equivalent to 10mm depth (density = 100 kg/m3)

		// NOTES: things to add or improve
		// 1. lake model: explicit representation of lake water storage, sunlight through lake
		//    with different purity, turbulent mixing of surface laker water, snow on frozen lake, etc.
		// 2. shallow snow wihtout a layer: melting energy
		// 3. urban model to be added.
		// 4. irrigation
		/// <summary>
		/// 刘永和添加的单位线选项，即hillslope汇流，坡面汇流
		/// 0为不计算汇流，相当于纯陆面模型
		/// 1为无单位延迟汇流，2为Nash单位线汇流，3为一半出流
		/// </summary>
		public static int OPT_UHG = 1;

		public static  int threads = 0;
		/// <summary>
		/// metric potential for wilting point (m)
		/// </summary>
		public static double PSIWLT = -150;
		
		
		public  int getYearLen(int year)
		{
			if (year % 400 == 0) {
				return 366;
			}
			if (year % 100 == 0) {
				return 365;
			}
			if (year % 4 == 0)
				return 366;
			
			return 365;
		}


		
		
		
		
		
		 
		

	}
}
