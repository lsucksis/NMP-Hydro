/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2018/1/14
 * Time: 19:16
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;
using System.Threading;

namespace NoahMP
{
	/// <summary>
	/// Description of NoahMP_Cell.
	/// </summary>
	public class GridCell
	{
		public	int ILOC;
		public int JLOC;
		public double LATITUDE;
		public double LONGITUDE;
		/// <summary>
		/// height
		/// </summary>
		public double HGT;
		public bool NeedCal = true;
		public int YEARLEN;
		public int JULIAN;
		/// <summary>
		/// cosine solar zenith angle [0-1]
		/// </summary>
		public double COSZ = 0;
		
		/// <summary>
		/// thickness of lowest layer
		/// </summary>
		public double DZ8W = 30;
		/// <summary>
		/// maximum no. of soil layers
		/// </summary>
		//public int NSOIL;
		/// <summary>
		/// layer-bottom depth from soil surf (m  <0)
		/// </summary>
		public FortDoubleArray ZSOIL;
		//		/// <summary>
		//		/// maximum no. of snow layers
		//		/// </summary>
		//		public int NSNOW;
		// // IN : Model configuration
		
		public double RSMIN {
			get {
				return REDPRM.RSTBL[VEGTYP - 1];
			}
		}
		public double HS {
			get {
				return REDPRM.HSTBL[VEGTYP - 1];
			}
		}
		public double RGL {
			get {
				return REDPRM.RGLTBL[VEGTYP - 1];
			}
		}
		//		public double SHDFAC {
		//			get {
		//				return REDPRM.SHDTBL[this.VEGTYP-1];
		//			}
		//		}
		/// <summary>
		/// yearly max vegetation fraction
		/// </summary>
		public double SHDMAX;
		public double SHDMIN;
		public double KDT {
			get {
				return NoahMP.REFKDT * this.DKSAT / NoahMP.REFDK;
			}
		}
		/// <summary>
		/// vegetation type
		/// </summary>
		public int VEGTYP;
		public int SOILTYP;
		public int NROOT {
			get {
				return REDPRM.NROTBL[this.VEGTYP - 1];
			}
		}
		//public int ISURBAN;
		/// <summary>
		/// ice (sea ice = 1,land ice =-1, neither =0)
		/// </summary>
		public int ICE;
		/// <summary>
		/// surface type 1->soil; 2->lake
		/// </summary>
		public int IST = 1;
		
		/// <summary>
		/// soil color type (1-lighest; 8-darkest)
		/// </summary>
		public int ISC;
		public FortDoubleArray SMCEQ;
		// // IN : Vegetation/Soil characteristics
		public int IZ0TLND;
		public double SFCTMP;
		// // IN : User options
		/// <summary>
		/// pressure (pa)
		/// </summary>
		public double SFCPRS;
		/// <summary>
		/// pressure at lowest model layer
		/// </summary>
		public double PSFC;
		public double UU;
		public double VV;
		/// <summary>
		/// mixing ratio (kg/kg)
		/// </summary>
		public double Q2;
		/// <summary>
		/// cloud water mixing ratio
		/// </summary>
		public double QC;
		/// <summary>
		/// downward shortwave radiation (w/m2)
		/// </summary>
		public double SOLDN;
		/// <summary>
		/// downward longwave radiation (w/m2)
		/// </summary>
		public double LWDN;
		/// <summary>
		/// precipitation rate (kg m-2 s-1)
		/// </summary>
		public double PRCP;
		/// <summary>
		/// bottom condition for soil temp. (k)
		/// </summary>
		public double TBOT = 273;
		/// <summary>
		/// Depth (m) of lower boundary soil temperature ( Assigned in REDPRM )
		/// </summary>
		public double ZBOT = -8;
		/// <summary>
		/// atmospheric co2 concentration (pa)
		/// </summary>
		public double CO2AIR {
			get {
				double CO2 = 395e-06;
				return CO2 * this.SFCPRS;
			}
		}
		/// <summary>
		/// atmospheric o2 concentration (pa)
		/// </summary>
		public double O2AIR {
			get {
				
				double O2 = 0.209;
				return O2 * this.SFCPRS;
			}
			
		}
		/// <summary>
		/// foliage nitrogen (%)
		/// </summary>
		public double FOLN;
		/// <summary>
		/// ice fraction at last timestep
		/// </summary>
		public FortDoubleArray FICEOLD;
		/// <summary>
		/// planetary boundary layer height
		/// </summary>
		public double PBLH = 0;
		/// <summary>
		/// reference height (m)
		/// </summary>
		public double ZLVL;
		/// <summary>
		/// snow albedo at last time step (CLASS type)
		/// </summary>
		public double ALBOLD;
		/// <summary>
		/// snow mass at last time step (mm)
		/// </summary>
		public double SNEQVO;
		
		/// <summary>
		/// snow/soil temperature [k]
		/// </summary>
		public FortDoubleArray STC;
		/// <summary>
		/// liquid soil moisture [m3/m3]，在PHASECHANGE函数中修改
		/// </summary>
		public FortDoubleArray SH2O = null;
		/// <summary>
		/// soil moisture (ice + liq.) [m3/m3]
		/// </summary>
		public FortDoubleArray SMC = null;
		/// <summary>
		/// canopy air temperature (k)
		/// </summary>
		public double TAH = -9999;
		// 263.901062;
		/// <summary>
		/// canopy air vapor pressure (pa)
		/// </summary>
		public double EAH = -9999;
		// 930.417236;
		// 202.858627;
		/// <summary>
		/// fraction of canopy that is wet
		/// </summary>
		public double FWET;
		/// <summary>
		/// intercepted liquid water (mm)
		/// </summary>
		public double CANLIQ;
		public double CANICE;
		/// <summary>
		/// vegetation temperature (k)
		/// </summary>
		public double TV = -9999;
		//268.471619;// 268.872559;// 265.540253;
		/// <summary>
		/// ground temperature (k)
		/// </summary>
		public double TG = -9999;
		//268.872559;// 265.540253;
		
		/// <summary>
		/// !mixing ratio at lowest model layer
		/// </summary>
		public double QSFC = 100000;
		/// <summary>
		/// snow at ground srf (mm/s) [+]
		/// </summary>
		public double QSNOW = 0e-6;
		
		/// <summary>
		/// 最上一层雪层的序号:0-2。原fortran版为默认的三层雪序号为 -2,-1,0, 四层土壤序号为1,2,3,4
		/// </summary>
		public int ISNOW;
		/// <summary>
		/// depth of snow/soil layer-bottom
		/// </summary>
		public FortDoubleArray ZSNSO;
		/// <summary>
		/// snow height [m]
		/// </summary>
		public double SNOWH = 0;
		/// <summary>
		/// snow water eqv. [mm]
		/// </summary>
		public double SNEQV;
		/// <summary>
		/// snow layer ice [mm]
		/// </summary>
		public FortDoubleArray SNICE = new FortDoubleArray(1 - Driver.NSnow, 0);
		//new double[Driver.NSnow];
		/// <summary>
		/// snow layer liquid water [mm]
		/// </summary>
		public FortDoubleArray SNLIQ = new FortDoubleArray(1 - Driver.NSnow, 0);
		/// <summary>
		/// the depth to water table [m]
		/// </summary>
		public double ZWT = 2.5;
		/// <summary>
		/// water storage in aquifer [mm]
		/// </summary>
		public double WA = 4900;
		/// <summary>
		/// water storage in aquifer
		/// </summary>
		public double WT = 4900;
		/// <summary>
		/// lake water storage (can be neg.) (mm)
		/// </summary>
		public double WSLAKE;
		/// <summary>
		/// leaf mass [g/m2]
		/// </summary>
		public double LFMASS = 50;
		/// <summary>
		/// mass of fine roots [g/m2]
		/// </summary>
		public double RTMASS = 500;
		/// <summary>
		/// stem mass [g/m2]
		/// </summary>
		public double STMASS = 50;
		/// <summary>
		/// mass of wood (incl. woody roots) [g/m2]
		/// </summary>
		public double WOOD = 500;
		/// <summary>
		/// stable carbon in deep soil [g/m2]
		/// </summary>
		public double STBLCP = 1000;
		/// <summary>
		/// short-lived carbon, shallow soil [g/m2]
		/// </summary>
		public double FASTCP = 1000;
		/// <summary>
		/// leaf area index [-]
		/// </summary>
		public double LAI;
		/// <summary>
		/// stem area index [-]
		/// </summary>
		public double SAI = 0.1;
		/// <summary>
		/// momentum drag coefficient
		/// </summary>
		public double CM = 0.1;
		/// <summary>
		/// sensible heat exchange coefficient
		/// </summary>
		public double CH = 0.1;
		/// <summary>
		/// non-dimensional snow age
		/// </summary>
		public double TAUSS;
		/// <summary>
		/// soil water content between bottom of the soil and water table [m3/m3]
		/// </summary>
		public double SMCWTD;
		/// <summary>
		/// recharge to or from the water table when deep [m]
		/// </summary>
		public double DEEPRECH;
		/// <summary>
		/// recharge to or from the water table when shallow [m] (diagnostic)
		/// </summary>
		public double RECH;
		/// <summary>
		/// total absorbed solar radiation (w/m2)
		/// </summary>
		public double FSA;
		/// <summary>
		/// total reflected solar radiation (w/m2)
		/// </summary>
		public double FSR;
		/// <summary>
		/// total net LW rad (w/m2)  [+ to atm]
		/// </summary>
		public double FIRA;
		/// <summary>
		/// total sensible heat (w/m2) [+ to atm]
		/// </summary>
		public double FSH;
		/// <summary>
		/// /ground heat flux (w/m2)   [+ to soil]
		/// </summary>
		public double SSOIL;
		/// <summary>
		/// canopy evap heat (w/m2) [+ to atm]
		/// </summary>
		public double FCEV;
		/// <summary>
		/// ground evap heat (w/m2) [+ to atm]
		/// </summary>
		public double FGEV;
		/// <summary>
		/// transpiration heat (w/m2) [+ to atm]
		/// </summary>
		public double FCTR;
		/// <summary>
		/// evaporation of intercepted water (mm/s)
		/// </summary>
		public double ECAN;
		/// <summary>
		/// transpiration rate (mm/s)
		/// </summary>
		public double ETRAN;
		/// <summary>
		/// soil surface evaporation rate (mm/s]
		/// </summary>
		public double EDIR;
		/// <summary>
		/// surface radiative temperature (k)
		/// </summary>
		public double TRAD = 273;
		/// <summary>
		/// ground surface temp. [k]
		/// </summary>
		public double TGB = 273;
		/// <summary>
		/// ground surface temp. [k]
		/// </summary>
		public double TGV = 273;
		/// <summary>
		/// 2 m height air temperature (k)
		/// </summary>
		public double T2MV = -9999;
		/// <summary>
		/// 2 m height air temperature (k)
		/// </summary>
		public double T2MB = -9999;
		public double Q2V;
		/// <summary>
		/// bare ground heat conductance
		/// </summary>
		public double Q2B;
		/// <summary>
		/// surface runoff [mm/s]
		/// </summary>
		public double RUNSRF = 0;
		/// <summary>
		/// baseflow (sturation excess) [mm/s]
		/// </summary>
		public double RUNSUB = 0;
		public double AccRunSrf;
		public double AccRunSub;
		/// <summary>
		/// 坡面汇流，即进入河道汇流的入流量
		/// </summary>
		public double Q_HillSlope;
		/// <summary>
		/// par absorbed per unit LAI (w/m2)
		/// </summary>
		public double APAR;
		/// <summary>
		/// total leaf photosyn (umolco2/m2/s) [+]
		/// </summary>
		public double PSN;
		/// <summary>
		/// solar rad absorbed by veg. (w/m2)
		/// </summary>
		public double SAV;
		/// <summary>
		/// solar rad absorbed by ground (w/m2)
		/// </summary>
		public double SAG;
		/// <summary>
		/// snow cover fraction on the ground (-)
		/// </summary>
		public double FSNO;
		/// <summary>
		/// net ecosys exchange (g/m2/s CO2)
		/// </summary>
		public double NEE;
		/// <summary>
		/// net instantaneous assimilation [g/m2/s C]
		/// </summary>
		public double GPP;
		/// <summary>
		/// net primary productivity [g/m2/s C]
		/// </summary>
		public double NPP;
		/// <summary>
		/// 
		/// </summary>
		public double FVEG;
		public double FVEGMP;
		/// <summary>
		/// 反照率
		/// </summary>
		public double ALBEDO;
		/// <summary>
		/// melting water out of snow bottom [mm/s]
		/// </summary>
		public double QSNBOT;
		
		public double PONDING;
		public double PONDING1;
		public double PONDING2;
		/// <summary>
		/// sunlit leaf stomatal resistance (s/m)
		/// </summary>
		public double RSSUN;
		/// <summary>
		/// shaded leaf stomatal resistance (s/m)
		/// </summary>
		public double RSSHA;
		/// <summary>
		/// between canopy gap fraction for beam (-)
		/// </summary>
		public double BGAP;
		/// <summary>
		/// within canopy gap fraction for beam (-)
		/// </summary>
		public double WGAP;
		/// <summary>
		/// sensible heat exchange coefficient over vegetated fraction
		/// </summary>
		public double CHV;
		/// <summary>
		/// sensible heat exchange coefficient over bare-ground
		/// </summary>
		public double CHB;
		public double EMISSI;
		/// <summary>
		/// ground sen. heat [w/m2]   [+ to atm]
		/// </summary>
		public double SHG;
		/// <summary>
		/// canopy sen. heat [w/m2]   [+ to atm]
		/// </summary>
		public double SHC;
		/// <summary>
		/// sensible heat [w/m2]     [+ to atm]
		/// </summary>
		public double SHB;
		/// <summary>
		/// evaporation heat flux (w/m2)  [+= to atm]
		/// </summary>
		public double EVG;
		/// <summary>
		/// latent heat flux (w/m2)   [+ to atm]
		/// </summary>
		public double EVB;
		/// <summary>
		/// ground heat flux [w/m2]  [+ to soil]
		/// </summary>
		public double GHV;
		/// <summary>
		/// ground heat flux [w/m2] [+ to soil]
		/// </summary>
		public double GHB;
		/// <summary>
		/// ground net LW rad. [w/m2] [+ to atm]
		/// </summary>
		public double IRG;
		/// <summary>
		/// canopy net LW rad. [w/m2] [+ to atm]
		/// </summary>
		public double IRC;
		/// <summary>
		/// net longwave rad. [w/m2] [+ to atm]
		/// </summary>
		public double IRB;
		/// <summary>
		/// transpiration heat [w/m2] [+ to atm]
		/// </summary>
		public double TR;
		/// <summary>
		/// canopy evap. heat [w/m2]  [+ to atm]
		/// </summary>
		public double EVC;
		/// <summary>
		/// leaf exchange coefficient
		/// </summary>
		public double CHLEAF;
		/// <summary>
		/// under canopy exchange coefficient
		/// </summary>
		public double CHUC;
		/// <summary>
		/// sensible heat exchange coefficient over vegetated fraction
		/// </summary>
		public double CHV2;
		/// <summary>
		/// sensible heat exchange coefficient over bare-ground
		/// </summary>
		public double CHB2;
		/// <summary>
		/// snow fraction in precipitation
		/// </summary>
		public double FPICE;
		//#ifdef WRF_HYDRO
		public double SFCHEADRT;
		
		public double sumPcp;
		public		double sumEvap;
		public		double sumRunoff;	
				
		public double QFX {
			get {
				return ECAN + EDIR + ETRAN;
			}
		}
		public double CANWAT {
			get {
				return CANLIQ + CANICE;
			}
		}
		public double LH {
			get {
				return FCEV + FGEV + FCTR;
			}
		}
		/// <summary>
		/// B parameter ( Assigned in REDPRM )
		/// </summary>
		public double BEXP {
			get {
				return NoahMP.BB[this.SOILTYP - 1];
			}
		}

		/// <summary>
		/// soil thermal diffusivity/conductivity coef ( Assigned in REDPRM )
		/// </summary>
		public double F1 {
			get {
				return NoahMP.F11[this.SOILTYP - 1];
			}
		}
		/// <summary>
		/// 饱和土壤含水率
		/// </summary>
		public double SMCMAX {
			get {
				return NoahMP.MAXSMC[this.SOILTYP - 1];
			}
		}
		/// <summary>
		/// reference soil moisture (field capacity) (volumetric) ( Assigned in REDPRM )
		/// </summary>
		public double SMCREF {
			get {
				return NoahMP.REFSMC[this.SOILTYP - 1];
			}
		}
		public double FRZX {
			get {
				return NoahMP.FRZK * (SMCMAX / SMCREF) * (0.412 / 0.468);
			}
		}
		//
		/// <summary>
		/// saturated soil matric potential ( Assigned in REDPRM )
		/// </summary>
		public double PSISAT {
			get {
				return NoahMP.SATPSI[this.SOILTYP - 1];
			}
		}
		/// <summary>
		/// saturated soil hydraulic conductivity ( Assigned in REDPRM )
		/// </summary>
		public double DKSAT {
			get {
				return NoahMP.SATDK[this.SOILTYP - 1];
			}
		}
		/// <summary>
		/// saturated soil hydraulic diffusivity ( Assigned in REDPRM )
		/// </summary>
		public double DWSAT {
			get {
				return NoahMP.SATDW[this.SOILTYP - 1];
			}
		}
		/// <summary>
		/// wilting point soil moisture (volumetric) ( Assigned in REDPRM )
		/// </summary>
		public double SMCWLT {
			get {
				return NoahMP.WLTSMC[this.SOILTYP - 1];
			}
		}
		/// <summary>
		/// soil quartz content ( Assigned in REDPRM )
		/// </summary>
		public double QUARTZ {
			get {
				return NoahMP.QTZ[this.SOILTYP - 1];
			}
		}
		public int SLOPETYP = 2;
		/// <summary>
		/// 影响地下水的一个重要参数
		/// </summary>
		public double SLOPE {
			get {
				//return 0.1;
				return REDPRM.SLOPE_DATA[this.SLOPETYP-1];
			}
		}
		
		

		public double gamma(double z)
		{
			double a = z + 1 / (12 * z - 0.1 / z);
			return Math.Sqrt(2 * Math.PI / z) * Math.Pow(1 / Math.E * a, z);
		}
		public double NashUnit(double t)
		{
			double N = 3.5;
			double K = 3.0;
			return 1.0 / K / gamma(N) * Math.Pow(t / K, N - 1) * Math.Exp(-t / K);
		}

		
		/// <summary>
		/// 用于记录最近多个时段的Runoff，主要是服务于后续汇流的单位线卷积计算
		/// 0号值记录着当前时段的runoff，1号值记录着前一时段的，...
		/// </summary>
		public double[] HistoricalRunoff;
		/// <summary>
		/// 刘永和自己加的，用于调试
		/// </summary>
		double preValue = -1e20;
		/// <summary>
		/// 刘永和自己加的，用于调试
		/// </summary>
		int ncount=0;
		public void CheckModel()
		{
//			if(preValue>0.01 && this.FVEG<1e-15){
//				throw new Exception("FVEG become very zero:"+this.FVEG);
//			}
//			preValue=this.FVEG;
		}				
		

		public void	SurfaceFlux()
		{
			double COSZ = Driver.CALC_DECLIN(Driver.time0, this.LATITUDE, this.LONGITUDE);
			int YEARLEN = NoahMP.getYearLen(Driver.time0.Year);
//			if(TG<100)
//				throw new Exception("");
//			if(this.SHDMAX>1 || this.SHDMAX<0.01 || this.SHDMIN>1)
//				throw new Exception("");
			NoahMP.NoahMP_SFLX(this, -1, -1, this.LATITUDE, YEARLEN, -1, COSZ,
				Driver.DT, Driver.DX, DZ8W, Driver.NSoil, ZSOIL, Driver.NSnow,
				FVEG, SHDMAX, VEGTYP, Driver.vegparams.ISURBAN, ICE, IST,
				ISC, SMCEQ,
				IZ0TLND, SFCTMP, SFCPRS, PSFC, UU, VV, Q2,
				QC, SOLDN, LWDN, PRCP, TBOT, CO2AIR,
				O2AIR, FOLN, FICEOLD, PBLH, ZLVL,
				ref ALBOLD,ref SNEQVO,
				STC, SH2O, SMC,ref TAH,ref EAH,ref FWET,
				ref CANLIQ, ref CANICE, ref TV, ref TG, ref QSFC,ref QSNOW,
				ref ISNOW, ZSNSO, ref SNOWH,ref SNEQV, SNICE, SNLIQ,
				ref ZWT, ref WA, ref WT, ref WSLAKE, ref LFMASS, ref RTMASS,
				ref STMASS, ref WOOD, ref STBLCP, ref FASTCP, ref LAI, ref SAI,
				ref CM, ref CH, ref TAUSS,
				ref SMCWTD, ref DEEPRECH, ref RECH,
				out FSA, out FSR, out FIRA, out FSH, out SSOIL, out FCEV,
				out FGEV, out FCTR, out ECAN, out ETRAN, out EDIR, out TRAD,
				out TGB, out TGV, out T2MV, out T2MB, out Q2V, out Q2B,
				out RUNSRF, out RUNSUB, out APAR, out PSN, out SAV, out SAG,
				out FSNO, out NEE, out GPP, out NPP, out FVEGMP, out ALBEDO,
				out QSNBOT, out PONDING, out PONDING1, out PONDING2, out RSSUN, out RSSHA,
				out BGAP, out WGAP, out CHV, out CHB, out EMISSI,
				out SHG, out SHC, out SHB, out EVG, out EVB, out GHV,
				out GHB, out IRG, out IRC, out IRB, out TR, out EVC,
				out  CHLEAF, out CHUC, out CHV2, out CHB2, out FPICE,
				ref  SFCHEADRT);
//			Console.WriteLine(ALBEDO);
			
		}
		
		
		
		public GridCell()
		{
		}
	}
}
