/*
 * Created by SharpDevelop.
 * User: Yonghe Liu, Henan Polytechnic University
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
	/// Description of NoahMP.
	/// </summary>
	public class NoahMP3:NoahMP
	{		
		public static  void ATM(double SFCPRS, double SFCTMP, double Q2, double PRCP, double SOLDN, double COSZ,
			out double THAIR, out double QAIR, out double EAIR, out double RHOAIR,
			out double QPRECC, out double QPRECL, double[] SOLAD, double[] SOLAI,
			out double SWDOWN)
		{
			//PAIR   //atm bottom level pressure (pa)
			
			double PAIR = SFCPRS;                   // atm bottom level pressure (pa)
			THAIR = SFCTMP * Math.Pow(SFCPRS / PAIR, RAIR / CPAIR);

			//       QAIR   = Q2 / (1.0+Q2)           // mixing ratio to specific humidity [kg/kg]
			QAIR = Q2;                       // In WRF, driver converts to specific humidity

			EAIR = QAIR * SFCPRS / (0.622 + 0.378 * QAIR);
			RHOAIR = (SFCPRS - 0.378 * EAIR) / (RAIR * SFCTMP);

			//this.PRCP=PRCP=0;
			QPRECC = 0.10 * PRCP;          // should be from the atmospheric model
			QPRECL = 0.90 * PRCP;          // should be from the atmospheric model

			if (COSZ <= 0)
				SWDOWN = 0;
			else
				SWDOWN = SOLDN;
			
			SOLAD[0] = SWDOWN * 0.7 * 0.5;     // direct  vis
			SOLAD[1] = SWDOWN * 0.7 * 0.5;     // direct  nir
			SOLAI[0] = SWDOWN * 0.3 * 0.5;     // diffuse vis
			SOLAI[1] = SWDOWN * 0.3 * 0.5;     // diffuse nir

		}
		
		/// <summary>
		/// 
		/// </summary>
		/// <param name="cell"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="LAT"></param>
		/// <param name="YEARLEN"></param>
		/// <param name="JULIAN"></param>
		/// <param name="COSZ"></param>
		/// <param name="DT"></param>
		/// <param name="DX"></param>
		/// <param name="DZ8W"></param>
		/// <param name="NSOIL"></param>
		/// <param name="ZSOIL"></param>
		/// <param name="NSNOW"></param>
		/// <param name="SHDFAC"></param>
		/// <param name="SHDMAX"></param>
		/// <param name="VEGTYP"></param>
		/// <param name="ISURBAN"></param>
		/// <param name="ICE"></param>
		/// <param name="IST"></param>
		/// <param name="ISC"></param>
		/// <param name="SMCEQ"></param>
		/// <param name="IZ0TLND"></param>
		/// <param name="SFCTMP"></param>
		/// <param name="SFCPRS"></param>
		/// <param name="PSFC"></param>
		/// <param name="UU"></param>
		/// <param name="VV"></param>
		/// <param name="Q2"></param>
		/// <param name="QC"></param>
		/// <param name="SOLDN"></param>
		/// <param name="LWDN"></param>
		/// <param name="PRCP"></param>
		/// <param name="TBOT"></param>
		/// <param name="CO2AIR"></param>
		/// <param name="O2AIR"></param>
		/// <param name="FOLN">预留的氮限制作用</param>
		/// <param name="FICEOLD"></param>
		/// <param name="PBLH"></param>
		/// <param name="ZLVL"></param>
		/// <param name="ALBOLD"></param>
		/// <param name="SNEQVO"></param>
		/// <param name="STC"></param>
		/// <param name="SH2O"></param>
		/// <param name="SMC"></param>
		/// <param name="TAH"></param>
		/// <param name="EAH"></param>
		/// <param name="FWET"></param>
		/// <param name="CANLIQ"></param>
		/// <param name="CANICE"></param>
		/// <param name="TV"></param>
		/// <param name="TG"></param>
		/// <param name="QSFC"></param>
		/// <param name="QSNOW"></param>
		/// <param name="ISNOW"></param>
		/// <param name="ZSNSO"></param>
		/// <param name="SNOWH"></param>
		/// <param name="SNEQV"></param>
		/// <param name="SNICE"></param>
		/// <param name="SNLIQ"></param>
		/// <param name="ZWT"></param>
		/// <param name="WA"></param>
		/// <param name="WT"></param>
		/// <param name="WSLAKE"></param>
		/// <param name="LFMASS"></param>
		/// <param name="RTMASS"></param>
		/// <param name="STMASS"></param>
		/// <param name="WOOD"></param>
		/// <param name="STBLCP"></param>
		/// <param name="FASTCP"></param>
		/// <param name="LAI"></param>
		/// <param name="SAI"></param>
		/// <param name="CM"></param>
		/// <param name="CH"></param>
		/// <param name="TAUSS"></param>
		/// <param name="SMCWTD"></param>
		/// <param name="DEEPRECH"></param>
		/// <param name="RECH"></param>
		/// <param name="FSA"></param>
		/// <param name="FSR"></param>
		/// <param name="FIRA"></param>
		/// <param name="FSH"></param>
		/// <param name="SSOIL"></param>
		/// <param name="FCEV"></param>
		/// <param name="FGEV"></param>
		/// <param name="FCTR"></param>
		/// <param name="ECAN"></param>
		/// <param name="ETRAN"></param>
		/// <param name="EDIR"></param>
		/// <param name="TRAD"></param>
		/// <param name="TGB"></param>
		/// <param name="TGV"></param>
		/// <param name="T2MV"></param>
		/// <param name="T2MB"></param>
		/// <param name="Q2V"></param>
		/// <param name="Q2B"></param>
		/// <param name="RUNSRF"></param>
		/// <param name="RUNSUB"></param>
		/// <param name="APAR"></param>
		/// <param name="PSN"></param>
		/// <param name="SAV"></param>
		/// <param name="SAG"></param>
		/// <param name="FSNO"></param>
		/// <param name="NEE"></param>
		/// <param name="GPP"></param>
		/// <param name="NPP"></param>
		/// <param name="FVEG"></param>
		/// <param name="ALBEDO"></param>
		/// <param name="QSNBOT"></param>
		/// <param name="PONDING"></param>
		/// <param name="PONDING1"></param>
		/// <param name="PONDING2"></param>
		/// <param name="RSSUN"></param>
		/// <param name="RSSHA"></param>
		/// <param name="BGAP"></param>
		/// <param name="WGAP"></param>
		/// <param name="CHV"></param>
		/// <param name="CHB"></param>
		/// <param name="EMISSI"></param>
		/// <param name="SHG"></param>
		/// <param name="SHC"></param>
		/// <param name="SHB"></param>
		/// <param name="EVG"></param>
		/// <param name="EVB"></param>
		/// <param name="GHV"></param>
		/// <param name="GHB"></param>
		/// <param name="IRG"></param>
		/// <param name="IRC"></param>
		/// <param name="IRB"></param>
		/// <param name="TR"></param>
		/// <param name="EVC"></param>
		/// <param name="CHLEAF"></param>
		/// <param name="CHUC"></param>
		/// <param name="CHV2"></param>
		/// <param name="CHB2"></param>
		/// <param name="FPICE"></param>
		/// <param name="SFCHEADRT"></param>
		public static void NoahMP_SFLX(GridCell cell,
			int ILOC, int  JLOC, double LAT, int YEARLEN, double JULIAN, double COSZ, // ! IN : Time/Space-related
			double DT, double DX, double DZ8W, int NSOIL, FortDoubleArray ZSOIL, int NSNOW, // ! IN : Model configuration
			double SHDFAC, double SHDMAX, int VEGTYP, int ISURBAN, int ICE, int IST, // ! IN : Vegetation/Soil characteristics
			int ISC, FortDoubleArray SMCEQ,                                         // ! IN : Vegetation/Soil characteristics
			int IZ0TLND,                                                   // ! IN : User options
			double SFCTMP, double SFCPRS, double PSFC, double UU, double VV, double Q2, // ! IN : Forcing
			double QC, double SOLDN, double LWDN, double PRCP, double TBOT, double CO2AIR, // ! IN : Forcing
			double  O2AIR, double FOLN, FortDoubleArray  FICEOLD, double PBLH, double ZLVL,           // ! IN : Forcing
			ref double ALBOLD, ref double SNEQVO,                                         // ! IN/OUT :
			FortDoubleArray STC, FortDoubleArray SH2O, FortDoubleArray SMC, ref double TAH, ref double EAH, ref double FWET, // ! IN/OUT :
			ref double CANLIQ, ref double CANICE, ref double TV, ref double TG, ref double QSFC, ref double QSNOW, // ! IN/OUT :
			ref int ISNOW, FortDoubleArray ZSNSO, ref double SNOWH, ref double SNEQV, FortDoubleArray SNICE, FortDoubleArray SNLIQ, // ! IN/OUT :
			ref double ZWT, ref double WA, ref double WT, ref double WSLAKE, ref double LFMASS, ref double RTMASS, // ! IN/OUT :
			ref double STMASS, ref double WOOD, ref double STBLCP, ref double FASTCP, ref double LAI, ref double SAI, // ! IN/OUT :
			ref double CM, ref double CH, ref double TAUSS,                               // ! IN/OUT :
			ref double  SMCWTD, ref double DEEPRECH, ref double RECH,                               // ! IN/OUT :
			out double FSA, out double FSR, out double FIRA, out double FSH, out double SSOIL, out double FCEV, // ! OUT :
			out double FGEV, out double FCTR, out double ECAN, out double ETRAN, out double EDIR, out double TRAD, // ! OUT :
			out double  TGB, out double   TGV, out double   T2MV, out double   T2MB, out double   Q2V, out double Q2B, // ! OUT :
			out double   RUNSRF, out double   RUNSUB, out double APAR, out double PSN, out double SAV, out double SAG, // ! OUT :
			out double   FSNO, out double   NEE, out double   GPP, out double NPP, out double FVEG, out double ALBEDO, // ! OUT :
			out double   QSNBOT, out double   PONDING, out double   PONDING1, out double   PONDING2, out double   RSSUN, out double   RSSHA, // ! OUT :
			out double   BGAP, out double   WGAP, out double   CHV, out double   CHB, out double   EMISSI,           // ! OUT :
			out double   SHG, out double   SHC, out double   SHB, out double   EVG, out double   EVB, out double   GHV, // ! OUT :
			out double    GHB, out double   IRG, out double   IRC, out double   IRB, out double   TR, out double   EVC, // ! OUT :
			out double    CHLEAF, out double   CHUC, out double CHV2, out double CHB2, out double FPICE,               //
			ref double SFCHEADRT   // ! IN/OUT :
		)
		{
			
			QC = -1e36;                     // test dummy value
			PBLH = -1E36;                     // test dummy value // PBL height
			NEE = 0;
			GPP = 0;
			NPP = 0;
			ALBEDO = 0;
			//this.DZ8W1D = 2 * ZLVL;                          // thickness of atmospheric layers
			//COSZ = Driver.CALC_DECLIN(Driver.time0, LATITUDE, LONGITUDE);
//			if(double.IsNaN(TGB) || double.IsNaN(TG))
//						throw new Exception("");
			// only if DVEG == 2.
			
			// --------------------------------------------------------------------------------------------------
			// Re-process atmospheric forcing
			double RHOAIR = 0; //kg/m3 
			double QAIR = 0;//specific humidity (kg/kg) (q2/(1+q2))
			double QPRECC = 0;
			double QPRECL = 0;
			double SWDOWN = 0;
			double THAIR = 0;
			double EAIR = 0;
			double[] SOLAD = new double[2];
			double[] SOLAI = new double[2];
			ATM(SFCPRS, SFCTMP, Q2, PRCP, SOLDN, COSZ, out THAIR,
				out QAIR, out EAIR, out RHOAIR, out QPRECC, out QPRECL, SOLAD, SOLAI,
				out SWDOWN);
			// snow/soil layer thickness (m)
			
			FortDoubleArray DZSNSO = new FortDoubleArray(1 - NSNOW, NSOIL);
			for (int IZ = ISNOW + 1; IZ <= Driver.NSoil; IZ++) {
				if (IZ == ISNOW + 1)
					DZSNSO[IZ] = -ZSNSO[IZ];
				else
					DZSNSO[IZ] = ZSNSO[IZ - 1] - ZSNSO[IZ];
			}
			
			
			// root-zone temperature
			double TROOT = 0;
			//int NROOT =	REDPRM.NROTBL[this.VEGTYP-1];
			for (int IZ = 1; IZ <= cell.NROOT; IZ++) {
				//STC--snow/soil temperature  DZSNSO--snow/soil layer thickness [m]
				TROOT += STC[IZ] * DZSNSO[IZ] / (-ZSOIL[cell.NROOT]);
			}
			
			// Total water storage for water balance check
			if (IST == 1) {
				double BEG_WB = CANLIQ + CANICE + SNEQV + WA;
				for (int IZ = 1; IZ <= Driver.NSoil; IZ++)
					BEG_WB += SMC[IZ] * DZSNSO[IZ] * 1000;
			}

//			if (double.IsNaN(SNOWH))
//				throw new Exception();
			// Vegetation phenology
			double HTOP = 0;
			double ELAI = 0; //LAI在雪埋后的校正
			double ESAI = 0;
			double IGS = 0;
			//对DVEG为动态植被时，PHENOLOGY不对LAI、SAI做修改，LAI请查看后面的CARBON函数
			PHENOLOGY(cell, VEGTYP, Driver.vegparams.ISURBAN, SNOWH, TV, cell.LATITUDE, YEARLEN, JULIAN, ref LAI, ref SAI, TROOT,
				out HTOP, out ELAI, out ESAI, out IGS);
//			if (double.IsNaN(ELAI))
//				throw new Exception();
			FVEG = 0;
			if (NoahMP.DVEG == 1) {
				FVEG = SHDFAC;//静态不更新
				if (FVEG <= 0.01)
					FVEG = 0.01;
			} else if (NoahMP.DVEG == 2 || NoahMP.DVEG == 3) {
				FVEG = 1 - Math.Exp(-0.52 * (LAI + SAI));   //动态植被的关键行
				if (FVEG <= 0.01)
					FVEG = 0.01;
			} else if (NoahMP.DVEG == 4 || NoahMP.DVEG == 5) {
				FVEG = SHDMAX;
				if (FVEG <= 0.01)
					FVEG = 0.01;
			} else {
				throw new Exception("Wrong DVEG value " + DVEG);
			}
			
			if (VEGTYP == Driver.vegparams.ISURBAN || VEGTYP == Driver.vegparams.ISBARREN)
				FVEG = 0;
			if (Math.Abs(ELAI + ESAI) < 1e-20)
				FVEG = 0; 

			/// <summary>
			/// phase change index [1-melt; 2-freeze]
			/// </summary>
			FortIntArray IMELT = new FortIntArray(1 - Driver.NSnow, Driver.NSoil);
			FortDoubleArray SNICEV = new FortDoubleArray(1 - Driver.NSnow, 0);
			//  new double[Driver.NSnow];
			//partial volume ice of snow [m3/m3]
			FortDoubleArray SNLIQV = new FortDoubleArray(1 - Driver.NSnow, 0);
			//  new double[Driver.NSnow];
			//partial volume liq of snow [m3/m3]
			FortDoubleArray EPORE = new FortDoubleArray(1 - Driver.NSnow, 0);
			double T2M = 0;
			double QMELT = 0;
			double TAUX = 0;
			double TAUY = 0;
			//soil water evaporation factor (0 - 1)
			FortDoubleArray BTRANI = new FortDoubleArray(1, Driver.NSoil);
			double BTRAN = 0;
			double TS = 0;
			double LATHEAV = 0;
			//latent heat vap./sublimation (j/kg)
			double LATHEAG = 0;
			//latent heat vap./sublimation (j/kg)
			bool FROZEN_GROUND = false;
			// used to define latent heat pathway
			bool FROZEN_CANOPY = false;
			double FSRV = 0;
			//
			double FSRG = 0;
			//double ERRWAT = 0;
			//water error [kg m{-2}]
			double Q1 = 0;
			double Q2E = 0;
			/// <summary>
			/// soil ice content (m3/m3)
			/// </summary>
			FortDoubleArray SICE = new FortDoubleArray(1, Driver.NSoil);
			Energy(cell, ICE, VEGTYP, IST, ISC, Driver.NSnow, Driver.NSoil,  //in
				ISNOW, cell.NROOT, Driver.DT, RHOAIR, SFCPRS, QAIR,  //in
				SFCTMP, THAIR, LWDN, UU, VV, ZLVL,  //in
				CO2AIR, O2AIR, SOLAD, SOLAI, COSZ, IGS,  //in
				EAIR, HTOP, TBOT, ZBOT, ZSNSO, ZSOIL,  //in
				ELAI, ESAI, NoahMP.CSOIL, FWET, FOLN,          //in
				FVEG,                                          //in
				QSNOW, DZSNSO, cell.LATITUDE, CANLIQ, CANICE, ILOC, JLOC,  //in
				IMELT, SNICEV, SNLIQV, EPORE, out T2M, out FSNO,  //out
				out SAV, out SAG, QMELT, out FSA, out FSR, out TAUX,  //out
				out TAUY, out FIRA, out FSH, out FCEV, out FGEV, out FCTR,  //out
				out TRAD, out PSN, out APAR, out SSOIL, BTRANI, out BTRAN,  //out
				out PONDING, out TS, out LATHEAV, out LATHEAG, out FROZEN_CANOPY, out FROZEN_GROUND,                          //ref
				ref TV, ref TG, STC, ref SNOWH, ref EAH, ref TAH,  //inref
				ref SNEQVO, ref SNEQV, SH2O, SMC, SNICE, SNLIQ,  //inref
				ref ALBOLD, ref CM, ref CH, Driver.DX, ref DZ8W, ref Q2,  //inref
				ref TAUSS,                                          //inref
				ref QC, ref  PBLH, ref QSFC, ref PSFC, Driver.vegparams.ISURBAN, IZ0TLND,  //in
				out T2MV, out T2MB, out FSRV,
				out FSRG, out RSSUN, out RSSHA, out BGAP, out WGAP, out TGV, out TGB,
				out Q1, out Q2V, out Q2B, out Q2E, out CHV, out CHB,  //out
				out EMISSI, out SHG, out SHC, out SHB, out EVG, out EVB, out GHV, out GHB, out IRG, out IRC, out IRB,
				out TR, out EVC, out CHLEAF, out CHUC, out CHV2, out CHB2);  //out
			
			//估算土壤冰含量
			for (int i = 1; i <= Driver.NSoil; i++) {
				SICE[i] = Math.Max(0.0, SMC[i] - SH2O[i]);
			}
			SNEQVO = SNEQV;

			
			double QVAP = Math.Max(FGEV / LATHEAG, 0);       // positive part of fgev; Barlage change to ground v3.6
			double QDEW = Math.Abs(Math.Min(FGEV / LATHEAG, 0));  // negative part of fgev
			EDIR = QVAP - QDEW;

			// Compute water budgets (water storages, ET components, and runoff)

			//WATER函数中的run选项已完全通过，除了5
			double CMC = 0;
			/// <summary>
			/// groundwater recharge [mm/s]
			/// </summary>
			double QIN = 0;
		
			/// <summary>
			/// groundwater discharge [mm/s]
			/// </summary>
			double QDIS = 0;
			WATER(cell, VEGTYP, Driver.NSnow, Driver.NSoil, IMELT, Driver.DT, UU,  //in
				VV, FCEV, FCTR, QPRECC, QPRECL, ELAI,  //in
				ESAI, SFCTMP, QVAP, QDEW, ZSOIL, BTRANI,  //in
				FICEOLD, PONDING, TG, IST, FVEG, ILOC, JLOC, SMCEQ,  //in
				LATHEAV, LATHEAG, FROZEN_CANOPY, FROZEN_GROUND,                         //in  MB
				ref ISNOW, ref CANLIQ, ref CANICE, ref TV, ref SNOWH, ref SNEQV,  //inout
				SNICE, SNLIQ, STC, ZSNSO, SH2O, SMC,  //inout
				SICE, ref ZWT, ref WA, ref WT, DZSNSO, ref WSLAKE,  //inout
				ref	SMCWTD, ref DEEPRECH, ref RECH,  //inout
				out	CMC, out ECAN, out ETRAN, out FWET, out RUNSRF, out RUNSUB,  //out
				out	QIN, out QDIS, out QSNOW, out PONDING1, out PONDING2,
				Driver.vegparams.ISURBAN, out QSNBOT, out FPICE,
			      //#ifdef WRF_HYDRO
				SFCHEADRT);  //out
			//下面这句为刘永和所加
			CalHillSlopeFlow(cell, RUNSRF + RUNSUB);
			
			//Console.WriteLine("SH2O="+SH2O[0]);
			//WATER已基本通过连续多个时次的模拟验证
			// Compute carbon budgets (carbon storages and co2  bvoc fluxes)
			if (NoahMP.DVEG == 2 || NoahMP.DVEG == 5) {
				//2021年时已验证过一次，结果一致
				double AUTORS = 0;
				//net ecosystem respiration (g/m2/s C)
				double HETERS = 0;
				//organic respiration (g/m2/s C)
				double TOTSC = 0;
				//  //total soil carbon (g/m2)
				double TOTLB = 0;
				//total living carbon (g/m2)
				CARBON(Driver.NSnow, Driver.NSoil, VEGTYP, cell.NROOT, Driver.DT, ZSOIL,  //in
					DZSNSO, STC, SMC, TV, TG, PSN,  //in
					FOLN, cell.SMCMAX, BTRAN, APAR, FVEG, IGS,  //in
					TROOT, IST, cell.LATITUDE, ILOC, JLOC, Driver.vegparams.ISURBAN,  //in
					ref	LFMASS, ref RTMASS, ref STMASS, ref WOOD, ref STBLCP, ref FASTCP,  //inout
					out GPP, out NPP, out NEE, out AUTORS, out HETERS, out TOTSC,  //out
					out TOTLB, ref LAI, ref SAI);                   //out
				//Console.WriteLine("LAI" + LAI + " SAI" + SAI);
			}

			// water and energy balance check

			//     ERROR (SWDOWN ,FSA    ,FSR    ,FIRA   ,FSH    ,FCEV   ,  //in
			//                 FGEV   ,FCTR   ,SSOIL  ,BEG_WB ,CANLIQ ,CANICE ,  //in
			//                 SNEQV  ,WA     ,SMC    ,DZSNSO ,PRCP   ,ECAN   ,  //in
			//                 ETRAN  ,EDIR   ,RUNSRF ,RUNSUB ,DT     ,NSOIL  ,  //in
			//                 NSNOW  ,IST    ,ERRWAT ,ILOC   , JLOC  ,FVEG   ,
			//                 SAV    ,SAG    ,FSRV   ,FSRG   ,ZWT  );   //in ( Except ERRWAT, which is out )
			// urban - jref
			double QFX = ETRAN + ECAN + EDIR;
			if (VEGTYP == Driver.vegparams.ISURBAN) {
				QSFC = (QFX / RHOAIR * CH) + QAIR;
				Q2B = QSFC;
			}

			if (SNOWH <= 1E-6 || SNEQV <= 1E-3) {
				SNOWH = 0.0;
				SNEQV = 0.0;
			}

			if (SWDOWN > 1e-20)
				ALBEDO = FSR / SWDOWN;
			else
				ALBEDO = -1e20;

			
		}
		public static void WATER(GridCell cell, int VEGTYP, int NSNOW, int NSOIL, FortIntArray IMELT, double DT, double UU, double VV,   //in
			double FCEV, double FCTR, double QPRECC, double QPRECL, double ELAI, double ESAI, double SFCTMP, double QVAP,   //in
			double QDEW, FortDoubleArray ZSOIL, FortDoubleArray BTRANI, FortDoubleArray FICEOLD, double PONDING, double TG,  //in
			int IST, double FVEG, int ILOC, int JLOC, FortDoubleArray  SMCEQ, double LATHEAV, double  LATHEAG,  //in
			bool  frozen_canopy, bool FROZEN_GROUND,                         //in  MB
			ref int  ISNOW, ref double CANLIQ, ref double CANICE, ref double TV, ref double SNOWH, ref double SNEQV,    //inout
			FortDoubleArray SNICE, FortDoubleArray SNLIQ, FortDoubleArray STC, FortDoubleArray ZSNSO, FortDoubleArray SH2O, FortDoubleArray SMC,    //inout
			FortDoubleArray SICE, ref double ZWT, ref double WA, ref double WT, FortDoubleArray DZSNSO, ref double WSLAKE,    //inout
			ref double SMCWTD, ref double DEEPRECH, ref double RECH,   //inout
			out  double  CMC, out double ECAN, out double ETRAN, out double FWET, out double RUNSRF, out double RUNSUB,    //out
			out double QIN, out double QDIS, out double QSNOW, out double PONDING1, out double PONDING2,
			int ISURBAN, out double QSNBOT, out double FPICE,	double sfcheadrt)
		{
			double QRAIN = 0;   //rain at ground srf (mm) [+]
			double SNOWHIN = 0;// //snow depth increasing rate (m/s)
			FortDoubleArray ETRANI = new FortDoubleArray(1, Driver.NSoil);  //transpiration rate (mm/s) [+]
			FortDoubleArray WCND = new FortDoubleArray(1, Driver.NSoil);   //hydraulic conductivity (m/s)
			double QDRAIN = 0;//soil-bottom free drainage [mm/s]
			double SNOFLOW = 0; //glacier flow [mm/s]
			double FCRMAX = 0; //maximum of FCR (-)
			double WSLMAX = 5000;      //maximum lake water storage (mm)

			// initialize
			
			QIN = 0;
			QDIS = 0;
			QSNBOT = 0;
			SNOFLOW = 0;
			
			double QINSUR = 0;  //water input on soil surface [m/s]
			// canopy-intercepted snowfall/rainfall, drips, and throughfall

			//下面的CANWATER函数没有严格测试过ISNOW<3的情况
			//敏感参数：TG,FCTR
			CANWATER(VEGTYP, DT, SFCTMP, UU, VV,
				FCEV, FCTR, QPRECC, QPRECL, ELAI,
				ESAI, IST, TG, FVEG, ILOC, JLOC,
				frozen_canopy,
				ref CANLIQ, ref CANICE, ref TV,
				out CMC, out ECAN, out ETRAN, out QRAIN, out QSNOW,
				out SNOWHIN, out FWET, out FPICE);
			

			// sublimation, frost, evaporation, and dew

			double QSNSUB = 0; //snow surface sublimation rate [mm/s]
			if (SNEQV > 0) {
				QSNSUB = Math.Min(QVAP, SNEQV / DT);
			}
			double QSEVA = QVAP - QSNSUB; //QSEVA: soil surface evap rate [mm/s]

			double QSNFRO = 0;
			if (SNEQV > 0) {
				QSNFRO = QDEW;
			}
			if (double.IsNaN(QDEW))
				throw new Exception();
			double QSDEW = QDEW - QSNFRO;

			//函数调用会更新SH2O的顶层
			SNOWWATER(Driver.NSnow, Driver.NSoil, IMELT, DT, ZSOIL,
				SFCTMP, SNOWHIN, QSNOW, QSNFRO, QSNSUB,
				QRAIN, FICEOLD, ILOC, JLOC,
				ref ISNOW, ref SNOWH, ref SNEQV, SNICE, SNLIQ,
				SH2O, SICE, STC, ZSNSO, DZSNSO,
				out QSNBOT, out SNOFLOW, out PONDING1, out PONDING2);

			
			if (FROZEN_GROUND) {
				SICE[1] += (QSDEW - QSEVA) * DT / (DZSNSO[1] * 1000);
				QSDEW = 0.0;
				QSEVA = 0.0;
				if (SICE[1] < 0) {					
					SH2O[1] += SICE[1];
					SICE[1] = 0;
				}
			}
			// convert units (mm/s -> m/s)

			//PONDING: melting water from snow when there is no layer
			QINSUR = (PONDING + PONDING1 + PONDING2) / DT * 0.001;   //    QINSUR = PONDING/DT * 0.001
//			if(double.IsNaN(QSNBOT))
//				throw new Exception("");
			if (ISNOW == 0)
				QINSUR += (QSNBOT + QSDEW + QRAIN) * 0.001;
			else
				QINSUR += (QSNBOT + QSDEW) * 0.001;
			

			QSEVA *= 0.001;
			
			for (int iz = 1; iz <= cell.NROOT; iz++) {
				ETRANI[iz] = ETRAN * BTRANI[iz] * 0.001;
			}

			QINSUR += cell.SFCHEADRT / DT * 0.001;   //sfcheadrt units (m)

			

			RUNSUB = 0;
			RUNSRF = 0;
			if (IST == 2) { // lake
				RUNSRF = 0;
				if (WSLAKE >= WSLMAX)
					RUNSRF = QINSUR * 1000;   //mm/s
				WSLAKE += (QINSUR - QSEVA) * 1000 * DT - RUNSRF * DT;   //mm
				//RUNSUB = 0;
			} else {			
				//如果不是湖，则运行下面的				
//				Console.WriteLine("SH2O[4]=" + SH2O[4]);
				SOILWATER(cell, NSOIL, NSNOW, DT, ZSOIL, DZSNSO,
					QINSUR, QSEVA, ETRANI, SICE, ILOC, JLOC,
					SH2O, SMC, ref ZWT, ISURBAN, VEGTYP,
					ref SMCWTD, ref DEEPRECH,
					out RUNSRF, out QDRAIN, out RUNSUB, WCND, out FCRMAX);

				if (NoahMP.OPT_RUN == 1) {					
					//已调试过，应该是对的了
					GROUNDWATER(cell, NSNOW, NSOIL, DT, SICE, ZSOIL,
						STC, WCND, FCRMAX, ILOC, JLOC,
						SH2O, ref ZWT, ref WA, ref WT, out QIN, out QDIS);
//					Console.Write("SICE=" + SICE[1].ToString("e") + "\t" + SICE[2].ToString("e") + "\t" + SICE[2].ToString("e") + "\t" + SICE[4].ToString("e") + "\t");
//					Console.Write("STC=" + STC[1].ToString("e") + "\t" + STC[2].ToString("e") + "\t" + STC[3].ToString("e") + "\t" + STC[4].ToString("e") + "\t");
//					Console.Write("WCND=" + WCND[1].ToString("e") + "\t" + WCND[2].ToString("e") + "\t" + WCND[3].ToString("e") + "\t" + WCND[4].ToString("e") + "\t");
//					Console.Write("FCRMAX=" + FCRMAX.ToString("e") + "\t");
//					Console.Write("SH2O=" + SH2O[1].ToString("e") + "\t" + SH2O[2].ToString("e") + "\t" + SH2O[3].ToString("e") + "\t" + SH2O[4].ToString("e") + "\t");
//					Console.WriteLine("ZWT=" + ZWT.ToString("e") + "\tWA=" + WA.ToString("e") + "\tWT=" + WT.ToString("e") + "\tQIN=" + QIN.ToString("e") + "\tQDIS=" + QDIS.ToString("e"));
					RUNSUB = QDIS;   //mm/s
				}

				if (NoahMP.OPT_RUN == 3 || NoahMP.OPT_RUN == 4) {
					RUNSUB += QDRAIN;   //mm/s
				}
				for (int iz = 1; iz <= Driver.NSoil; iz++) {
					SMC[iz] = SH2O[iz] + SICE[iz]; //此句是唯一修改SMC的地方
				}
				if (NoahMP.OPT_RUN == 5) {
					//还未验证，由于缺少参数文件，无法使用
					SHALLOWWATERTABLE(cell, NSNOW, NSOIL, ZSOIL, DT,
						DZSNSO, SMCEQ, ILOC, JLOC,
						SMC, ref ZWT, ref SMCWTD, out RECH, ref QDRAIN);   //inout

					SH2O[NSOIL] = SMC[NSOIL] - SICE[NSOIL];
					RUNSUB += QDRAIN;   //it really comes from void watertable, which is not called with the same frequency as the soil routines here
					WA = 0;				
				}
			}
			RUNSUB += SNOFLOW;   //mm/s			
		}
		public static void CANWATER(int VEGTYP, double DT, double SFCTMP, double UU, double VV,
			double FCEV, double FCTR, double QPRECC, double QPRECL, double ELAI,
			double ESAI, double IST, double TG, double FVEG, int ILOC, int JLOC,
			bool FROZEN_CANOPY, ref double CANLIQ, ref double CANICE, ref double TV,
			out double CMC, out double ECAN, out double ETRAN, out double QRAIN, out double QSNOW,
			out double SNOWHIN, out double FWET, out double FPICE)
		{
			
			// ------------------------ local variables ---------------------------
			double MAXSNO;   //canopy capacity for snow interception (mm)
			double MAXLIQ;   //canopy capacity for rain interception (mm)
			double FP = 0;   //fraction of the gridcell that receives precipitation
			double BDFALL;   //bulk density of snowfall (kg/m3)
			double QINTR = 0;   //interception rate for rain (mm/s)
			double QDRIPR = 0;   //drip rate for rain (mm/s)
			double QTHROR = 0;   //throughfall for rain (mm/s)
			double QINTS = 0;   //interception (loading) rate for snowfall (mm/s)
			double QDRIPS = 0;   //drip (unloading) rate for intercepted snow (mm/s)
			double QTHROS = 0;   //throughfall of snowfall (mm/s)
			double QEVAC;   //evaporation rate (mm/s)
			double QDEWC;   //dew rate (mm/s)
			double QFROC;   //frost rate (mm/s)
			double QSUBC;   //sublimation rate (mm/s)
			double FT;   //temperature factor for unloading rate
			double FV;   //wind factor for unloading rate
			double QMELTC;   //melting rate of canopy snow (mm/s)
			double QFRZC;   //refreezing rate of canopy liquid water (mm/s)
			double RAIN = 0;   //rainfall (mm/s)
			double SNOW = 0;   //snowfall (mm/s)
			double CANMAS = 0;   //total canopy mass (kg/m2)

			
			QRAIN = 0.0;
			QSNOW = 0.0;
			SNOWHIN = 0.0;
			//ECAN = 0.0;
			
			// partition precipitation into rain and snow.
			// Jordan (1991)
			
			FPICE = 0;
			if (NoahMP.OPT_SNF == 1) {
				if (SFCTMP > NoahMP.TFRZ + 2.5) {
					FPICE = 0;
				} else if (SFCTMP <= NoahMP.TFRZ + 0.5) {
					FPICE = 1.0;
				} else if (SFCTMP <= NoahMP.TFRZ + 2) {
					FPICE = 1 - (-54.632 + 0.2 * SFCTMP);
				} else {
					FPICE = 0.6;
				}
				
			}

			if (NoahMP.OPT_SNF == 2) {
				if (SFCTMP >= NoahMP.TFRZ + 2.2)
					FPICE = 0;
				else
					FPICE = 1.0;
				
			}

			if (NoahMP.OPT_SNF == 3) {
				if (SFCTMP >= NoahMP.TFRZ)
					FPICE = 0;
				else
					FPICE = 1.0;
			}

			;   // Hedstrom NR and JW Pomeroy (1998), Hydrol. Processes, 12, 1611-1625
			;   // fresh snow density

			BDFALL = Math.Min(120, 67.92 + 51.25 * Math.Exp((SFCTMP - NoahMP.TFRZ) / 2.59));   // Barlage: change to MIN in v3.6

			RAIN = (QPRECC + QPRECL) * (1 - FPICE);
			SNOW = (QPRECC + QPRECL) * FPICE;

			;   // fractional area that receives precipitation (see, Niu et al. 2005)

			if (QPRECC + QPRECL > 0) {
				FP = (QPRECC + QPRECL) / (10 * QPRECC + QPRECL);
			}
//			Console.WriteLine("Total:"+(QPRECC + QPRECL));
			
			// --------------------------- liquid water ------------------------------
			//double MAXLIQ;   //canopy capacity for rain interception (mm)
			MAXLIQ = Driver.vegparams.CH2OP[VEGTYP - 1] * (ELAI + ESAI);

			// average interception and throughfall

			if ((ELAI + ESAI) > 0) {
				QINTR = FVEG * RAIN * FP;   // interception capability
				QINTR = Math.Min(QINTR, (MAXLIQ - CANLIQ) / DT * (1 - Math.Exp(-RAIN * DT / MAXLIQ)));
				QINTR = Math.Max(QINTR, 0);
				QDRIPR = FVEG * RAIN - QINTR;
				QTHROR = (1 - FVEG) * RAIN;
			} else {
				QINTR = 0;
				QDRIPR = 0;
				QTHROR = RAIN;
			}
			;   // evaporation, transpiration, and dew

			double HVAP = NoahMP.HVAP;
			double HSUB = NoahMP.HSUB;
			
			if (!FROZEN_CANOPY) {                // Barlage: change to frozen_canopy
				ETRAN = Math.Max(FCTR / HVAP, 0);
				QEVAC = Math.Max(FCEV / HVAP, 0);
				QDEWC = Math.Abs(Math.Min(FCEV / HVAP, 0));
				QSUBC = 0;
				QFROC = 0;
			} else {
				ETRAN = Math.Max(FCTR / HSUB, 0);
				QEVAC = 0;
				QDEWC = 0;
				QSUBC = Math.Max(FCEV / HSUB, 0);
				QFROC = Math.Abs(Math.Min(FCEV / HSUB, 0));
			}
			// canopy water balance. for convenience allow dew to bring CANLIQ above
			// maxh2o or else would have to re-adjust drip

			QEVAC = Math.Min(CANLIQ / DT, QEVAC);
			CANLIQ = Math.Max(0, CANLIQ + (QINTR + QDEWC - QEVAC) * DT);
			if (CANLIQ <= 1E-06)
				CANLIQ = 0.0;

			// --------------------------- canopy ice ------------------------------
			// for canopy ice

			MAXSNO = 6.6 * (0.27 + 46 / BDFALL) * (ELAI + ESAI);
			if ((ELAI + ESAI) > 0) {
				QINTS = FVEG * SNOW * FP;
				QINTS = Math.Min(QINTS, (MAXSNO - CANICE) / DT * (1 - Math.Exp(-SNOW * DT / MAXSNO)));
				QINTS = Math.Max(QINTS, 0);
				FT = Math.Max(0.0, (TV - 270.15) / 1.87E5);
				FV = Math.Sqrt(UU * UU + VV * VV) / 1.56E5;
				QDRIPS = Math.Max(0, CANICE) * (FV + FT);
				QTHROS = (1.0 - FVEG) * SNOW + (FVEG * SNOW - QINTS);
			} else {
				QINTS = 0;
				QDRIPS = 0;
				QTHROS = SNOW;
			}

			QSUBC = Math.Min(CANICE / DT, QSUBC);
			CANICE = Math.Max(0, CANICE + (QINTS - QDRIPS) * DT + (QFROC - QSUBC) * DT);
			if (CANICE <= 1E-6)
				CANICE = 0;
			;   // wetted fraction of canopy

			if (CANICE > 0)
				FWET = Math.Max(0, CANICE) / Math.Max(MAXSNO, 1E-06);
			else
				FWET = Math.Max(0, CANLIQ) / Math.Max(MAXLIQ, 1E-06);
			
			FWET = Math.Pow(Math.Min(FWET, 1), 0.667);

			;   // phase change

			QMELTC = 0;
			QFRZC = 0;

			double TFRZ = NoahMP.TFRZ;
			if (CANICE > 1E-6 && TV > TFRZ) {
				QMELTC = Math.Min(CANICE / DT, (TV - TFRZ) * NoahMP.CICE * CANICE / NoahMP.DENICE / (DT * NoahMP.HFUS));
				CANICE = Math.Max(0, CANICE - QMELTC * DT);
				CANLIQ = Math.Max(0, CANLIQ + QMELTC * DT);
				TV = FWET * TFRZ + (1 - FWET) * TV;
			}

			if (CANLIQ > 1E-6 && TV < TFRZ) {
				QFRZC = Math.Min(CANLIQ / DT, (TFRZ - TV) * NoahMP.CWAT * CANLIQ / NoahMP.DENH2O / (DT * NoahMP.HFUS));
				CANLIQ = Math.Max(0, CANLIQ - QFRZC * DT);
				CANICE = Math.Max(0, CANICE + QFRZC * DT);
				TV = FWET * TFRZ + (1 - FWET) * TV;
			}
			
			// total canopy water
			CMC = CANLIQ + CANICE;

			// total canopy evaporation
			ECAN = QEVAC + QSUBC - QDEWC - QFROC;

			// rain or snow on the ground
			QRAIN = QDRIPR + QTHROR;
			QSNOW = QDRIPS + QTHROS;
			SNOWHIN = QSNOW / BDFALL;

//			Console.WriteLine("QRAIN:"+QRAIN+" ECAN"+ECAN);
			
			if (IST == 2 && TG > TFRZ) {
				QSNOW = 0;
				SNOWHIN = 0;
			}
		}
		
		/// <summary>
		/// 
		/// </summary>
		/// <param name="VEGTYP"></param>
		/// <param name="ISURBAN"></param>
		/// <param name="SNOWH">雪深</param>
		/// <param name="TV"></param>
		/// <param name="LAT"></param>
		/// <param name="YEARLEN"></param>
		/// <param name="JULIAN"></param>
		/// <param name="LAI">叶面积指数</param>
		/// <param name="SAI">茎面积指数</param>
		/// <param name="TROOT"></param>
		/// <param name="HTOP"></param>
		/// <param name="ELAI"></param>
		/// <param name="ESAI"></param>
		/// <param name="IGS"></param>
		public static void PHENOLOGY(GridCell cell, int VEGTYP, int ISURBAN, double SNOWH, double TV, double LAT, int YEARLEN, double JULIAN,
			ref	double LAI, ref double SAI, double TROOT, out double HTOP, out double ELAI, out double ESAI, out double  IGS)
		{
			//如果不采用动态植被,则查表求LAI和SAI
			if (NoahMP.DVEG == 1 || NoahMP.DVEG == 3 || NoahMP.DVEG == 4) {
				
				double DAY = Driver.time0.DayOfYear - 1 + Driver.time0.Hour / 24.0;
				;// Northern Hemisphere
				if (cell.LATITUDE < 0) {  	// Southern Hemisphere.  DAY is shifted by 1/2 year.
					DAY = (DAY + (0.5 * 365)) % 366;
				}
				double T = 12.0 * DAY / 366;
				int IT1 = (int)(T + 0.5);
				int IT2 = IT1 + 1;
				double WT1 = (IT1 + 0.5) - T;
				double WT2 = 1 - WT1;
				if (IT1 < 1)
					IT1 = 12;
				if (IT2 > 12)
					IT2 = 1;

				LAI = WT1 * Driver.vegparams.LAIM[VEGTYP - 1, (int)IT1 - 1] + WT2 * Driver.vegparams.LAIM[VEGTYP - 1, (int)IT2 - 1];
				SAI = WT1 * Driver.vegparams.SAIM[VEGTYP - 1, (int)IT1 - 1] + WT2 * Driver.vegparams.SAIM[VEGTYP - 1, (int)IT2 - 1];
				
			}
			//下面的check屏蔽掉后，可以避免FVEG逐渐变小的bug。为了这个bug，另一个改动的地方在STOMATA函数中
			//if (SAI < 0.05)
			//	SAI = 0.0;                  // MB: SAI CHECK, change to 0.05 v3.6
			//if (LAI < 0.05 || SAI < Math.Abs(1e-20))
			//	LAI = 0.0;  // MB: LAI CHECK

			if ((VEGTYP == Driver.vegparams.ISWATER) || (VEGTYP == Driver.vegparams.ISBARREN) || (VEGTYP == Driver.vegparams.ISSNOW) || (VEGTYP == Driver.vegparams.ISURBAN)) {
				LAI = 0;
				SAI = 0;
			}
			//buried by snow
			double hvb = Driver.vegparams.HVB[VEGTYP - 1]; //植被底部高度m
			double hvt = Driver.vegparams.HVT[VEGTYP - 1]; //植被顶部高度m
			double DB = Math.Min(Math.Max(SNOWH - hvb, 0), hvt - hvb); //植被净高
			double FB = DB / Math.Max(1E-06, hvt - hvb);

			if (hvt > 0 & hvt <= 1.0) {   //MB: change to 1.0 and 0.2 to reflect
				double SNOWHC = hvt * Math.Exp(-SNOWH / 0.2);             //      changes to HVT in MPTABLE
//				if (SNOWHC < 1e-10)
//					SNOWHC = 1e-10;
				FB = Math.Min(SNOWH, SNOWHC) / SNOWHC;
			}

			ELAI = LAI * (1 - FB);//雪埋后的校正
			ESAI = SAI * (1 - FB);
			
			//下面的check屏蔽掉后，可以避免FVEG逐渐变小的bug,为了这个bug，另一个改动的地方在STOMATA函数中
			//当ESAI和ELAI都是0时-->会引起FVEG被设为0-->会导致VEGE_FLUX不执行-->PSN持续为0-->在CO2FLUX中进一步导致LAI为0，无法反转
			if (ESAI < 0.05)
				ESAI = 0.0;                   // MB: ESAI CHECK, change to 0.05 v3.6
			if (ELAI < 0.05 || ESAI < Math.Abs(1e-20))
				ELAI = 0.0;  // MB: LAI CHECK
//			Console.WriteLine("TV=" + TV);
			if (TV > Driver.vegparams.TMIN[VEGTYP - 1])
				IGS = 1;
			else
				IGS = 0;
			//Console.WriteLine("TV="+TV+"	IGS="+IGS);
			HTOP = hvt;
//			if(ELAI+ESAI<1e-10)
//				throw new Exception();
		}
		public static void Energy(GridCell cell, int ICE, int VEGTYP, int IST, int ISC, int NSNOW, int NSOIL,  //in
			int ISNOW, int NROOT, double DT, double RHOAIR, double SFCPRS, double QAIR,  //in
			double SFCTMP, double THAIR, double LWDN, double UU, double VV, double ZREF,  //in
			double CO2AIR, double O2AIR, double[] SOLAD, double[] SOLAI, double COSZ, double IGS,  //in
			double EAIR, double HTOP, double TBOT, double ZBOT, FortDoubleArray ZSNSO, FortDoubleArray ZSOIL,  //in
			double ELAI, double ESAI, double CSOIL, double FWET, double FOLN,          //in
			double FVEG,                                          //in
			double QSNOW, FortDoubleArray DZSNSO, double LAT, double CANLIQ, double CANICE, int ILOC, int JLOC,  //in
			FortIntArray  IMELT, FortDoubleArray SNICEV, FortDoubleArray SNLIQV, FortDoubleArray EPORE, out double T2M, out double FSNO,  //out
			out double  SAV, out double SAG, double QMELT, out double FSA, out double FSR, out double TAUX,  //out
			out double TAUY, out double FIRA, out double FSH, out double FCEV, out double FGEV, out double FCTR,  //out
			out double TRAD, out double PSN, out double APAR, out double SSOIL, FortDoubleArray BTRANI, out double BTRAN,  //out
			out double PONDING, out double TS, out double LATHEAV, out double LATHEAG, out bool frozen_canopy, out bool frozen_ground,     //out
			ref double TV, ref double TG, FortDoubleArray STC, ref double SNOWH, ref double EAH, ref double TAH,    //inref
			ref double SNEQVO, ref double SNEQV, FortDoubleArray SH2O, FortDoubleArray SMC, FortDoubleArray SNICE, FortDoubleArray SNLIQ,  //inref
			ref double ALBOLD, ref double CM, ref double CH, double DX, ref double DZ8W, ref double Q2, ref double TAUSS,
			ref double QC, ref double PBLH, ref double QSFC, ref double PSFC, int ISURBAN, int IZ0TLND,
			out double T2MV,   //in
			out double T2MB, out double FSRV, out double  FSRG, out double RSSUN, out double RSSHA, out double BGAP, out double WGAP, out double TGV, out double TGB,
			out double Q1, out double Q2V, out double Q2B, out double Q2E, out double CHV, out double CHB, out  double  EMISSI, out double SHG,
			out  double SHC, out double SHB, out double EVG, out double EVB, out double GHV, out double GHB,
			out double IRG, out double IRC, out double IRB, out double TR, out double EVC, out double CHLEAF, out double CHUC,
			out double CHV2, out double CHB2)
		{
			frozen_canopy = false;
			frozen_ground = false;
			
			T2MV = -9999;
			double TAUXV = 0;
			double TAUYV = 0;
			Q2V = 0;
			CHV2 = 0;
			CHLEAF = 0;
			CHUC = 0;
			IRB = 0;
			SHB = 0;
			Q2B = 0;
			CHB2 = 0;
			IRG = 0;
			IRC = 0;
			SHG = 0;
			EVG = 0;
			SHC = 0;
			EVC = 0;
			TR = 0;
			GHV = 0;
			// wind speed at reference height: ur >= 1
			double UR = Math.Max(Math.Sqrt(UU * UU + VV * VV), 1);
			// vegetated or non-vegetated

			double VAI = ELAI + ESAI;
			bool VEG = false;
			if (VAI > 0)
				VEG = true;

			// ground snow cover fraction [Niu and Yang, 2007, JGR]

			
			double Z0 = 0.01;// Bare-soil roughness length (m) (i.e., under the canopy)
			FSNO = 0;
			if (SNOWH > 0) {
				double BDSNO = SNEQV / SNOWH;
				double FMELT = Math.Pow(BDSNO / 100, NoahMP.M);
				FSNO = Math.Tanh(SNOWH / (2.5 * Z0 * FMELT));
			}

			// ground roughness length

			double Z0MG = 0;// z0 momentum, ground (m)
					
			if (IST == 2) {
				if (TG <= NoahMP.TFRZ) {
					Z0MG = 0.01 * (1.0 - FSNO) + FSNO * NoahMP.Z0SNO;
				} else {
					Z0MG = 0.01;
				}
			} else {
				Z0MG = Z0 * (1.0 - FSNO) + FSNO * NoahMP.Z0SNO;
			}
			
			// roughness length and displacement height

			double ZPDG = SNOWH;
			double Z0M = 0;
			double ZPD = 0;
			if (VEG) {
				Z0M = Driver.vegparams.Z0MVT[VEGTYP - 1];
				ZPD = 0.65 * HTOP;
				if (SNOWH > ZPD)
					ZPD = SNOWH;
			} else {
				Z0M = Z0MG;
				ZPD = ZPDG;
			}

			
			double ZLVL = Math.Max(ZPD, HTOP) + ZREF;
			if (ZPDG >= ZLVL)
				ZLVL = ZPDG + ZREF;
			//     UR   = UR*LOG(ZLVL/Z0M)/LOG(10./Z0M)       //input UR is at 10m

			// canopy wind absorption coeffcient

			double CWP = Driver.vegparams.CWPVT[VEGTYP - 1];
			// Thermal properties of soil, snow, lake, and frozen soil

			FortDoubleArray DF = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortDoubleArray HCPCT = new FortDoubleArray(-NSNOW + 1, NSOIL);
			FortDoubleArray FACT = new FortDoubleArray(-NSNOW + 1, NSOIL);
			ThermoProp(cell, NSOIL, NSNOW, ISNOW, IST, DZSNSO,
				DT, SNOWH, SNICE, SNLIQ, CSOIL,
				SMC, SH2O, TG, STC, UR,                            //不改变SH2O的值
				LAT, Z0M, ZLVL, VEGTYP, ISURBAN,
				DF, HCPCT, SNICEV, SNLIQV, EPORE,
				FACT);
			// Solar radiation: absorbed  reflected by the ground and canopy
			//ThermoProp首次输出已验证

			double FSUN = 0;			
			double LAISUN = 0;// sunlit leaf area index (m2/m2)			
			double LAISHA = 0;// shaded leaf area index (m2/m2)			
			double PSNSUN = 0;// sunlit photosynthesis (umolco2/m2/s)			
			double PSNSHA = 0;// shaded photosynthesis (umolco2/m2/s)			
			double PARSUN = 0;// PAR absorbed per sunlit LAI (w/m2)			
			double PARSHA = 0;// PAR absorbed per shaded LAI (w/m2)
			Radiation(VEGTYP, IST, ISC, ICE, Driver.NSoil,
				SNEQVO, SNEQV, DT, COSZ, SNOWH,
				TG, TV, FSNO, QSNOW, FWET,
				ELAI, ESAI, SMC, SOLAD, SOLAI,
				FVEG, ILOC, JLOC,
				ref ALBOLD, ref TAUSS,
				out FSUN, out LAISUN, out LAISHA, out PARSUN, out PARSHA, // //out
				out SAV, out SAG, out FSR, out FSA, out FSRV, //
				out FSRG, out BGAP, out WGAP);            //out
			// vegetation and ground emissivity
			

			double EMV = 1 - Math.Exp(-(ELAI + ESAI) / 1.0);
			double EMG = RAD_PARAMS.EG[IST - 1] * (1 - FSNO) + 1.0 * FSNO;
			if (ICE == 1)
				EMG = 0.98 * (1 - FSNO) + 1.0 * FSNO;			
				
			

			// soil moisture factor controlling stomatal resistance

			
			BTRAN = 0;
			//int NROOT = REDPRM.NROTBL[this.VEGTYP-1];
			double MPE = 1e-6;
			if (IST == 1) {
				double GX = 0;
				for (int IZ = 1; IZ <= NROOT; IZ++) {					
					if (NoahMP.OPT_BTR == 1) {                   // Noah
						GX = (SH2O[IZ] - cell.SMCWLT) / (cell.SMCREF - cell.SMCWLT);
					}
					if (NoahMP.OPT_BTR == 2) {                   // CLM
						double PSI = Math.Max(NoahMP.PSIWLT, -cell.PSISAT * Math.Pow(Math.Max(0.01, SH2O[IZ]) / cell.SMCMAX, -cell.BEXP));
						GX = (1 - PSI /NoahMP.PSIWLT) / (1 + cell.PSISAT / PSIWLT);
					}
					if (NoahMP.OPT_BTR == 3) {                   // SSiB
						double PSI = Math.Max(PSIWLT, -cell.PSISAT * Math.Pow(Math.Max(0.01, SH2O[IZ]) / cell.SMCMAX, -cell.BEXP));
						GX = 1 - Math.Exp(-5.8 * (Math.Log(PSIWLT / PSI)));
					}

					GX = Math.Min(1, Math.Max(0, GX));
					BTRANI[IZ] = Math.Max(MPE, DZSNSO[IZ] / (-ZSOIL[NROOT]) * GX);
					BTRAN += BTRANI[IZ];
				}
				BTRAN = Math.Max(MPE, BTRAN);
				for (int IZ = 1; IZ <= NROOT; IZ++)
					BTRANI[IZ] /= BTRAN;
			}
			//根据Science(2020)论文“Recent global decline of CO2 fertilization effects on vegetation photosynthesis”，GPP随CO2浓度的变化率beta为：每增加
            // 100ppm的CO2，Beta值减小0.7%/年。主要原因是由于(1)生物量变大后，氮磷肥料的限制造成的，(2)水分的限制。
			
			// soil surface resistance for ground evap.
			double RSURF = 0;
			double RHSUR = 0;
			double BEVAP = Math.Max(0.0, SH2O[1] / cell.SMCMAX);
			if (IST == 2) {
				RSURF = 1;          // avoid being divided by 0
				RHSUR = 1.0;
			} else {

				// RSURF based on Sakaguchi and Zeng, 2009
				// taking the "residual water content" to be the wilting point,
				// and correcting the exponent on the D term (typo in SZ09 ?)
				//L_RSURF = (-ZSOIL[1]) * ( exp ( (1.0 - Math.Min(1.0,SH2O[1]/SMCMAX)) ** 5 ) - 1.0 ) / ( 2.71828 - 1.0 );
				double L_RSURF = (-ZSOIL[1]) * (Math.Exp(Math.Pow(1.0 - Math.Min(1.0, SH2O[1] / cell.SMCMAX), 5)) - 1.0) / (2.71828 - 1.0);
				
				double D_RSURF = 2.2E-5 * cell.SMCMAX * cell.SMCMAX * Math.Pow(1.0 - cell.SMCWLT / cell.SMCMAX, 2.0 + 3.0 / cell.BEXP);
				RSURF = L_RSURF / D_RSURF;

				// Older RSURF computations:
				//    RSURF = FSNO * 1. + (1-FSNO)* Math.Exp(8.25-4.225*BEVAP) //Sellers (1992)
				//    RSURF = FSNO * 1. + (1-FSNO)* Math.Exp(8.25-6.0  *BEVAP) //adjusted to decrease RSURF for wet soil

				if (SH2O[1] < 0.01 && SNOWH <1e-8) {
					RSURF = 1E6;
				}
				double PSI = -cell.PSISAT * Math.Pow(Math.Max(0.01, SH2O[1]) / cell.SMCMAX, -cell.BEXP);
				RHSUR = FSNO + (1 - FSNO) * Math.Exp(PSI * NoahMP.GRAV / (NoahMP.RW * TG));
				//double rhsur=FSNO + (1 - FSNO) * Math.Exp(PSI * NoahMP.GRAV / (NoahMP.RW * 280)); //刘永和所加的调试语句
//				if (double.IsInfinity(RHSUR))
//					throw new Exception();
			}

			// urban - jref
			if (VEGTYP == ISURBAN && SNOWH <1e-8) {
				RSURF = 1E6;
			}
			// set psychrometric constant

			if (TV > NoahMP.TFRZ) {            // Barlage: add distinction between ground and
				LATHEAV = NoahMP.HVAP;                // vegetation in v3.6
				
				frozen_canopy = false;
			} else {
				LATHEAV = NoahMP.HSUB;
				frozen_canopy = true;
			}
			double GAMMAV = NoahMP.CPAIR * SFCPRS / (0.622 * LATHEAV);

			if (TG > NoahMP.TFRZ) {
				LATHEAG = NoahMP.HVAP;
				frozen_ground = false;
			} else {
				LATHEAG = NoahMP.HSUB;
				frozen_ground = true;
			}
			double GAMMAG = NoahMP.CPAIR * SFCPRS / (0.622 * LATHEAG);
			
			// Surface temperatures of the ground and canopy and energy fluxes
			TGV = 0;
			CHV = 0;
			RSSUN = 0;
			RSSHA = 0;

			double CMV = 0;
			if (VEG && FVEG > 0 && VEGTYP != Driver.vegparams.ISBARREN) {
				TGV = TG;
				CMV = CM;
				CHV = CH;
				VEGE_FLUX(cell, NSNOW, NSOIL, ISNOW, VEGTYP, VEG,  //in
					DT, SAV, SAG, LWDN, UR,  //in
					UU, VV, SFCTMP, THAIR, QAIR,  //in
					EAIR, RHOAIR, SNOWH, VAI, GAMMAV, GAMMAG,  //in
					FWET, LAISUN, LAISHA, CWP, DZSNSO,  //in
					HTOP, ZLVL, ZPD, Z0M, FVEG,  //in
					Z0MG, EMV, EMG, CANLIQ,  //in
					CANICE, STC, DF, out RSSUN, out RSSHA,
					RSURF, LATHEAV, LATHEAG, PARSUN, PARSHA, IGS,  //in
					FOLN, CO2AIR, O2AIR, BTRAN, SFCPRS,  //in
					RHSUR, ILOC, JLOC, Q2,  //in
					ref EAH, ref TAH, ref TV, ref TGV, ref CMV,  //inout   
					ref CHV, DX, DZ8W,                    //inout
					out TAUXV, out TAUYV, out IRG, out IRC, out SHG,  //out   
					out SHC, out EVG, out EVC, out TR, out GHV,  //out
					out T2MV, out PSNSUN, out PSNSHA,                    //out
				          //jref:start
					QC, PBLH, ref QSFC, PSFC, ISURBAN,  //in                      
					IZ0TLND, out Q2V, out CHV2, out CHLEAF, out CHUC);   
			} else {
//				Console.WriteLine("vege_flux no execution PSN="+(PSNSUN+PSNSHA));
				
				
			}
			
			TGB = TG;
			double CMB = CM;
			CHB = CH;
			double TAUXB = 0;
			double TAUYB = 0;
//			if (double.IsNaN(TGB) || double.IsNaN(SNOWH))
//				throw new Exception("");
			//此句之前已全部调试过，应该已无问题
			BARE_FLUX(NSNOW, NSOIL, ISNOW, DT, SAG,  //in
				LWDN, UR, UU, VV, SFCTMP,  //in
				THAIR, QAIR, EAIR, RHOAIR, SNOWH,  //in
				DZSNSO, ZLVL, ZPDG, Z0MG,           //in
				EMG, STC, DF, RSURF, LATHEAG,  //in
				GAMMAG, RHSUR, ILOC, JLOC, Q2,  //in
				ref TGB, ref CMB, ref CHB,                    //inout
				out TAUXB, out TAUYB, out IRB, out SHB, out EVB,  //out
				out GHB, out T2MB, DX, DZ8W, VEGTYP,  //out

				QC, PBLH, ref QSFC, PSFC, ISURBAN,  //in
				IZ0TLND, SFCPRS, out Q2B, out CHB2);                          //in
			//BARE_FLUX函数已测试通过
			
			//energy balance at vege canopy: SAV          =(IRC+SHC+EVC+TR)     *FVEG  at   FVEG
			//energy balance at vege ground: SAG*    FVEG =(IRG+SHG+EVG+GHV)    *FVEG  at   FVEG
			//energy balance at bare ground: SAG*(1-FVEG)=(IRB+SHB+EVB+GHB)*(1-FVEG) at 1-FVEG
			if (VEG && FVEG > 0) {
				TAUX = FVEG * TAUXV + (1.0 - FVEG) * TAUXB;
				TAUY = FVEG * TAUYV + (1.0 - FVEG) * TAUYB;
				FIRA = FVEG * IRG + (1.0 - FVEG) * IRB + IRC;
//				if (double.IsNaN(FIRA) || FIRA < 0)
//					throw new Exception("");
				FSH = FVEG * SHG + (1.0 - FVEG) * SHB + SHC;
				FGEV = FVEG * EVG + (1.0 - FVEG) * EVB;
				SSOIL = FVEG * GHV + (1.0 - FVEG) * GHB;
				FCEV = EVC;
				FCTR = TR;
				TG = FVEG * TGV + (1.0 - FVEG) * TGB;
//				if(double.IsNaN(TG))
//					throw new Exception("");
				T2M = FVEG * T2MV + (1.0 - FVEG) * T2MB;
				TS = FVEG * TV + (1.0 - FVEG) * TGB;
				CM = FVEG * CMV + (1.0 - FVEG) * CMB;      // better way to average?
				CH = FVEG * CHV + (1.0 - FVEG) * CHB;
				Q1 = FVEG * (EAH * 0.622 / (SFCPRS - 0.378 * EAH)) + (1.0 - FVEG) * QSFC;
				Q2E = FVEG * Q2V + (1.0 - FVEG) * Q2B;
			} else {
				TAUX = TAUXB;
				TAUY = TAUYB;
				FIRA = IRB;
				FSH = SHB;
				FGEV = EVB;
				SSOIL = GHB;
				TG = TGB;
//				if(double.IsNaN(TG))
//					throw new Exception("");
				T2M = T2MB;
				FCEV = 0;
				FCTR = 0;
				TS = TG;
				CM = CMB;
				CH = CHB;
				Q1 = QSFC;
				Q2E = Q2B;
				RSSUN = 0.0;
				RSSHA = 0.0;
				TGV = TGB;
				CHV = CHB;
			}

			double FIRE = LWDN + FIRA;
			if (TG < 100) {
				throw new Exception("");
			}

			if (FIRE <= 0) {
				//       WRITE(6,*) 'emitted longwave <0; skin T may be wrong due to inconsistent'
				//       WRITE(6,*) 'input of SHDFAC with LAI'
				//       WRITE(6,*) ILOC, JLOC, 'SHDFAC=',FVEG,'VAI=',VAI,'TV=',TV,'TG=',TG
				//       WRITE(6,*) 'LWDN=',LWDN,'FIRA=',FIRA,'SNOWH=',SNOWH
				//       wrf_error_fatal("STOP in Noah-MP")
//				throw new Exception("能量错误");
			}

			// Compute a net emissivity
			EMISSI = FVEG * (EMG * (1 - EMV) + EMV + EMV * (1 - EMV) * (1 - EMG)) + (1 - FVEG) * EMG;

			// When we're computing a TRAD, subtract from the emitted IR the
			// reflected portion of the incoming LWDN, so we're just
			// considering the IR originating in the canopy/ground system.

			//TRAD = ( ( FIRE - (1-EMISSI)*LWDN ) / (EMISSI*SB) ) ** 0.25;
			TRAD = Math.Pow((FIRE - (1 - EMISSI) * LWDN) / (EMISSI * NoahMP.SB), 0.25);
			// Old TRAD calculation not taking into account Emissivity:
			// TRAD = (FIRE/SB)**0.25
			//光合有效辐射 Photosynthetically active radiation
			APAR = PARSUN * LAISUN + PARSHA * LAISHA;
			PSN = PSNSUN * LAISUN + PSNSHA * LAISHA;
//			if(PSN>0)
//				throw new Exception();
			
			// 3L snow  4L soil temperatures
			TSNOSOI(ICE, NSOIL, NSNOW, ISNOW, IST,  //in
				TBOT, ZSNSO, SSOIL, DF, HCPCT,  //in
				ZBOT, SAG, DT, SNOWH, DZSNSO,  //in
				TG, ILOC, JLOC,                    //in
				STC);
//			for (int i = 0; i < NSNOW; i++) {
//				if (STC[i] < 200)
//					STC[i] = 273;
//			}
			//inout
			//上面的TSNOSOI基本测试正确，有较小的误差
			// adjusting snow surface temperature
			if (NoahMP.OPT_STC == 2) {
				if (SNOWH > 0.05 && TG > NoahMP.TFRZ) {
					TGV = NoahMP.TFRZ;
					TGB = NoahMP.TFRZ;
					if (VEG && FVEG > 0) {
						TG = FVEG * TGV + (1.0 - FVEG) * TGB;
						TS = FVEG * TV + (1.0 - FVEG) * TGB;
					} else {
						TG = TGB;
						TS = TGB;
					}
				}
			}

			PHASECHANGE(cell, NSNOW, NSOIL, ISNOW, DT, FACT, // //in
				DZSNSO, HCPCT, IST, ILOC, JLOC, //& //in
				STC, SNICE, SNLIQ, ref SNEQV, ref SNOWH, // !inout
				SMC, SH2O,                            // !inout
				out QMELT, IMELT, out PONDING);                   //  !out
			//PHASECHANGE测试通过
		}
		/// <summary>
		/// 计算水位深度平衡，通过土壤厚度，土层深度和体积含水量. 仅在RUN_OPT为2时被调用
		/// </summary>
		/// <param name="NSOIL"></param>
		/// <param name="NSNOW"></param>
		/// <param name="ZSOIL"></param>
		/// <param name="DZSNSO"></param>
		/// <param name="SH2O"></param>
		/// <returns></returns>
		public static double ZWTEQ(GridCell cell, int NSOIL, int NSNOW, FortDoubleArray ZSOIL, FortDoubleArray DZSNSO, FortDoubleArray SH2O)
		{
			
			// ----------------------------------------------------------------------
			// calculate equilibrium water table depth (Niu et al., 2005)
			// ----------------------------------------------------------------------
			//  IMPLICIT NONE
			// ----------------------------------------------------------------------
			// input//
			//  INTEGER,                         INTENT(IN)  NSOIL  //no. of soil layers
			//  INTEGER,                         INTENT(IN)  NSNOW  //maximum no. of snow layers
			//  double, DIMENSION(1:NSOIL),        INTENT(IN)  ZSOIL  //depth of soil layer-bottom [m]
			//  double, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  DZSNSO //snow/soil layer depth [m]
			//  double, DIMENSION(1:NSOIL),        INTENT(IN)  SH2O   //soil liquid water content [m3/m3]
			// output
			//  double,                           INTENT(OUT)  ZWT    //water table depth [m]
			
			//locals

			
			int K;
			; //do-loop index
			
			//double WD1; //water deficit from coarse (4-L) soil moisture profile
			//double WD2; //water deficit from fine (100-L) soil moisture profile
			//double DZFINE; //layer thickness of the 100-L soil layers to 6.0 m
			//double TEMP; //temporary variable			
			
			// ----------------------------------------------------------------------
			double WD1 = 0;
			for (K = 1; K <= NSOIL; K++) {
				WD1 += (cell.SMCMAX - SH2O[K]) * DZSNSO[K]; // [m]
			}
			
			int NFINE = 100; //no. of fine soil layers of 6m soil
			FortDoubleArray ZFINE = new FortDoubleArray(1, NFINE); //layer-bottom depth of the 100-L soil layers to 6.0 m
			double DZFINE = 3.0 * (-ZSOIL[NSOIL]) / NFINE;
			//do K =1,NFINE
			for (K = 1; K <= NFINE; K++) {
				ZFINE[K] = K * DZFINE;
			}

			double ZWT = -3 * ZSOIL[NSOIL] - 0.001;   // initial value [m]

			double WD2 = 0;
			for (K = 1; K <= NFINE; K++) {
				double TEMP = 1 + (ZWT - ZFINE[K]) / cell.PSISAT;
				WD2 += cell.SMCMAX * (1 - Math.Pow(TEMP, -1 / cell.BEXP)) * DZFINE;
				if (Math.Abs(WD2 - WD1) <= 0.01) {
					ZWT = ZFINE[K];
					break;
				}
			}
			return ZWT;
		}
		static void COMBO(ref double DZ, ref double  WLIQ, ref double  WICE, ref double T, double DZ2, double WLIQ2, double WICE2, double T2)
		{
			// input

			//    REAL, INTENT(IN)    :: DZ2   //nodal thickness of 2 elements being combined [m]
			//    REAL, INTENT(IN)    :: WLIQ2 //liquid water of element 2 [kg/m2]
			//    REAL, INTENT(IN)    :: WICE2 //ice of element 2 [kg/m2]
			//    REAL, INTENT(IN)    :: T2    //nodal temperature of element 2 [k]
			//    REAL, INTENT(INOUT) :: DZ    //nodal thickness of 1 elements being combined [m]
			//    REAL, INTENT(INOUT) :: WLIQ  //liquid water of element 1
			//    REAL, INTENT(INOUT) :: WICE  //ice of element 1 [kg/m2]
			//    REAL, INTENT(INOUT) :: T     //node temperature of element 1 [k]

			// local
//
			//    REAL                :: DZC   //total thickness of nodes 1 and 2 (DZC=DZ+DZ2).
			//    REAL                :: WLIQC //combined liquid water [kg/m2]
			//    REAL                :: WICEC //combined ice [kg/m2]
			//    REAL                :: TC    //combined node temperature [k]
			//    REAL                :: H     //enthalpy of element 1 [J/m2]
			//    REAL                :: H2    //enthalpy of element 2 [J/m2]
			//    REAL                :: HC    //temporary

			//-----------------------------------------------------------------------

			double DZC = DZ + DZ2;
			double WICEC = (WICE + WICE2);
			double WLIQC = (WLIQ + WLIQ2);
			double H = (NoahMP.CICE * WICE + NoahMP.CWAT * WLIQ) * (T - NoahMP.TFRZ) + NoahMP.HFUS * WLIQ;
			double H2 = (NoahMP.CICE * WICE2 + NoahMP.CWAT * WLIQ2) * (T2 - NoahMP.TFRZ) + NoahMP.HFUS * WLIQ2;
			double HC = H + H2;
			double TC;
			if (HC < 0) {
				TC = NoahMP.TFRZ + HC / (NoahMP.CICE * WICEC + NoahMP.CWAT * WLIQC);
			} else if (HC <= NoahMP.HFUS * WLIQC) {
				TC = NoahMP.TFRZ;
			} else {
				TC = NoahMP.TFRZ + (HC - NoahMP.HFUS * WLIQC) / (NoahMP.CICE * WICEC + NoahMP.CWAT * WLIQC);
			}


			DZ = DZC;
			WICE = WICEC;
			WLIQ = WLIQC;
			T = TC;

		}
		//END SUBROUTINE COMBO
		/// <summary>
		/// 核查过了，没问题,2022年11月16日
		/// </summary>
		/// <param name="NSNOW"></param>
		/// <param name="NSOIL"></param>
		/// <param name="ISNOW"></param>
		/// <param name="STC"></param>
		/// <param name="SNICE"></param>
		/// <param name="SNLIQ"></param>
		/// <param name="DZSNSO"></param>
		static void DIVIDE(int NSNOW, int NSOIL,                         // !in
			ref int ISNOW, FortDoubleArray STC, FortDoubleArray SNICE, FortDoubleArray SNLIQ, FortDoubleArray DZSNSO)  //inout
		{

			// input

			//    INTEGER, INTENT(IN)                            :: NSNOW //maximum no. of snow layers [ =3]
			//    INTEGER, INTENT(IN)                            :: NSOIL //no. of soil layers [ =4]

			// input and output

			//    INTEGER                        , INTENT(INOUT) :: ISNOW //actual no. of snow layers
			//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC   //snow layer temperature [k]
			//    REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE //snow layer ice [mm]
			//    REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ //snow layer liquid water [mm]
			//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO//snow layer depth [m]

			// local variables:

			//    INTEGER                                        :: J     //indices
			//    INTEGER                                        :: MSNO  //number of layer (top) to MSNO (bot)
			//    REAL                                           :: DRR   //thickness of the combined [m]
			FortDoubleArray DZ = new FortDoubleArray(1, NSNOW);    //snow layer thickness [m]
			FortDoubleArray SWICE = new FortDoubleArray(1, NSNOW); //partial volume of ice [m3/m3]
			FortDoubleArray SWLIQ = new FortDoubleArray(1, NSNOW); //partial volume of liquid water [m3/m3]
			FortDoubleArray TSNO = new FortDoubleArray(1, NSNOW);  //node temperature [k]
			//    REAL                                           :: ZWICE //temporary
			//    REAL                                           :: ZWLIQ //temporary
			//    REAL                                           :: PROPOR//temporary
			//    REAL                                           :: DTDZ  //temporary
			// ----------------------------------------------------------------------
			// 	DO J = 1,NSNOW
			for (int J = 1; J <= NSNOW; J++) {
				if (J <= Math.Abs(ISNOW)) {  //原句为J <= ABS(ISNOW)
					int index = J + ISNOW;
					DZ[J] = DZSNSO[index];  //原来：DZ的1，2层从DZSNCO的-1，0层取 现在：DZ的0，1层从DZSNCO的1，2层取
					SWICE[J] = SNICE[index];
					SWLIQ[J] = SNLIQ[index];
					TSNO[J] = STC[index];
				}
			}

			int	MSNO = Math.Abs(ISNOW);

			if (MSNO == 1) {  //IF (MSNO == 1) THEN/
				// Specify a new snow layer
				if (DZ[1] > 0.05) {
					MSNO = 2;
					DZ[1] /= 2;
					SWICE[1] /= 2;
					SWLIQ[1] /= 2;
					DZ[2] = DZ[1];
					SWICE[2] = SWICE[1];
					SWLIQ[2] = SWLIQ[1];
					TSNO[2] = TSNO[1];
				}
			}

			if (MSNO > 1) {
				if (DZ[1] > 0.05) {
					double	DRR = DZ[1] - 0.05;
					double	PROPOR = DRR / DZ[1];
					double ZWICE = PROPOR * SWICE[1];
					double ZWLIQ = PROPOR * SWLIQ[1];
					PROPOR = 0.05 / DZ[1];
					SWICE[1] = PROPOR * SWICE[1];
					SWLIQ[1] = PROPOR * SWLIQ[1];
					DZ[1] = 0.05;
					
					double dz = DZ[2];
					double swliq = SWLIQ[2];
					double swice = SWICE[2];
					double tsno = TSNO[2];
					COMBO(ref dz, ref swliq, ref swice, ref tsno, DRR, ZWLIQ, ZWICE, TSNO[1]);
					DZ[2] = dz;
					SWLIQ[2] = swliq;
					SWICE[2] = swice;
					TSNO[2] = tsno;
					// subdivide a new layer
					if (MSNO <= 2 && DZ[2] > 0.20) {  // MB: change limit
						//             if (MSNO <= 2 && DZ[1] > 0.10) {
						MSNO = 3;
						double DTDZ = (TSNO[1] - TSNO[2]) / ((DZ[1] + DZ[2]) / 2);
						DZ[2] /= 2;
						SWICE[2] /= 2;
						SWLIQ[2] /= 2;
						
						DZ[3] = DZ[2];
						SWICE[3] = SWICE[2];
						SWLIQ[3] = SWLIQ[2];
						TSNO[3] = TSNO[2] - DTDZ * DZ[2] / 2;
						if (TSNO[3] >= NoahMP.TFRZ) {
							TSNO[3] = TSNO[2];
						} else {
							TSNO[2] += DTDZ * DZ[2] / 2;
						}
					}
				}
			}

			if (MSNO > 2) { //IF (MSNO > 2) THEN
				if (DZ[2] > 0.2) {
					double DRR = DZ[2] - 0.2;
					double PROPOR = DRR / DZ[2];
					double ZWICE = PROPOR * SWICE[2];
					double ZWLIQ = PROPOR * SWLIQ[2];
					PROPOR = 0.2 / DZ[2];
					SWICE[2] *= PROPOR;
					SWLIQ[2] *= PROPOR;
					DZ[2] = 0.2;
					
					double dz = DZ[3];
					double swliq = SWLIQ[3];
					double swice = SWICE[3];
					double tsno = TSNO[3];
					COMBO(ref dz, ref swliq, ref swice, ref tsno, DRR, ZWLIQ, ZWICE, TSNO[2]);
					DZ[3] = dz;
					SWLIQ[3] = swliq;
					SWICE[3] = swice;
					TSNO[3] = tsno;
				}
			}

			ISNOW = -MSNO;
			//    DO J = ISNOW+1,0    
			for (int J = ISNOW + 1; J <= 0; J++) {
				int index = J - ISNOW;
				DZSNSO[J] = DZ[index];
				SNICE[J] = SWICE[index];
				SNLIQ[J] = SWLIQ[index];
				STC[J] = TSNO[index];
			}
		}
		/// <summary>
		/// 刘永和开发的函数，经过测试，这个平衡计算应该是不对的
		/// </summary>
		/// <param name="PDDUM"></param>
		/// <param name="DT"></param>
		/// <param name="ETRANI"></param>
		/// <param name="QSEVA"></param>
		/// <param name="SH2O"></param>
		/// <param name="RSAT"></param>
		/// <param name="DZSNSO"></param>
		/// <returns></returns>
		public static double TotalWater(double PDDUM, double DT, FortDoubleArray ETRANI, double QSEVA, FortDoubleArray SH2O, double RSAT, FortDoubleArray DZSNSO)
		{
			double sum0 = 0;
			double sum1 = 0;
			for (int i = 1; i <=NoahMP.NSOIL; i++) {				
				sum0 += SH2O[i] * DZSNSO[i];
				sum1 += ETRANI[i] * DT;
			}			
			return  RSAT * 1000 + PDDUM * DT - QSEVA * DT + sum0 - sum1;
		}
		/// <summary>
		/// 等同于NoahMP原版
		/// 注意，原代码中的两变量顺序是错误的，即(int ISURBAN, int VEGTYP)互换了，这里调整过来了
		/// 刘永和认为这个函数写得非常繁琐，可以很好地去简化
		/// </summary>
		/// <param name="NSOIL"></param>
		/// <param name="NSNOW"></param>
		/// <param name="DT"></param>
		/// <param name="ZSOIL"></param>
		/// <param name="DZSNSO"></param>
		/// <param name="QINSUR">water input on soil surface [mm/s]</param>
		/// <param name="QSEVA">土壤表面蒸发</param>
		/// <param name="ETRANI">土壤蒸发量mm/s</param>
		/// <param name="SICE"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="SH2O"></param>
		/// <param name="SMC"></param>
		/// <param name="ZWT"></param>
		/// <param name="ISURBAN"></param>
		/// <param name="VEGTYP"></param>
		/// <param name="SMCWTD">包气带土壤含水率</param>
		/// <param name="DEEPRECH"></param>
		/// <param name="RUNSRF"></param>
		/// <param name="QDRAIN">土壤底部排水</param>
		/// <param name="RUNSUB"></param>
		/// <param name="WCND">hydraulic conductivity (m/s)</param>
		/// <param name="FCRMAX"></param>
		public static void SOILWATER(GridCell cell, int NSOIL, int NSNOW, double DT, FortDoubleArray ZSOIL, FortDoubleArray DZSNSO,
			double QINSUR, double QSEVA, FortDoubleArray ETRANI, FortDoubleArray SICE, int ILOC, int JLOC,
			FortDoubleArray SH2O, FortDoubleArray SMC, ref double ZWT, int ISURBAN, int VEGTYP,
			ref	double SMCWTD, ref double DEEPRECH,
			out double RUNSRF, out double QDRAIN, out double RUNSUB, FortDoubleArray WCND, out double FCRMAX)
		{
			/*
			// input
			
  			int,                     INTENT(IN)  ILOC   //grid index
 			int,                     INTENT(IN)  JLOC   //grid index
  			int,                     INTENT(IN)  NSOIL  //no. of soil layers
  			int,                     INTENT(IN)  NSNOW  //maximum no. of snow layers
  			double,                        INTENT(IN)  DT     //time step (sec)
  			double, INTENT(IN)                         QINSUR //water input on soil surface [mm/s]
 			double, INTENT(IN)                         QSEVA  //evap from soil surface [mm/s]
  			double, DIMENSION(1:NSOIL),    INTENT(IN)  ZSOIL  //depth of soil layer-bottom [m]
  			double, DIMENSION(1:NSOIL),    INTENT(IN)  ETRANI //evapotranspiration from soil layers [mm/s]
  			double, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  DZSNSO //snow/soil layer depth [m]
  			double, DIMENSION(1:NSOIL), INTENT(IN)    SICE   //soil ice content [m3/m3]

  			int,                     INTENT(IN)  VEGTYP
  			int,                     INTENT(IN)  ISURBAN
			// input & output
  			double, DIMENSION(1:NSOIL), INTENT(INOUT)  SH2O   //soil liquid water content [m3/m3]
  			double, DIMENSION(1:NSOIL), INTENT(INOUT)  SMC    //total soil water content [m3/m3]
  			double, INTENT(INOUT)                      ZWT    //water table depth [m]
  			double,                     INTENT(INOUT)  SMCWTD //soil moisture between bottom of the soil and the water table [m3/m3]
  			double                    , INTENT(INOUT)  DEEPRECH

			// output
  			double, INTENT(OUT)                        QDRAIN //soil-bottom free drainage [mm/s]
  			double, INTENT(OUT)                        RUNSRF //surface runoff [mm/s]
  			double, INTENT(OUT)                        RUNSUB //subsurface runoff [mm/s]
  			double, INTENT(OUT)                        FCRMAX //maximum of FCR (-)
  			double, DIMENSION(1:NSOIL), INTENT(OUT)    WCND   //hydraulic conductivity (m/s)
			 */
			// local
			int K, IZ; //do-loop index
			int ITER; //iteration index
			double DTFINE; //fine time step (s)
			FortDoubleArray RHSTT = new FortDoubleArray(1, Driver.NSoil); //right-hand side term of the matrix DIMENSION(1:NSOIL)
			FortDoubleArray AI = new FortDoubleArray(1, Driver.NSoil); //left-hand side term, DIMENSION(1:NSOIL)
			FortDoubleArray BI = new FortDoubleArray(1, Driver.NSoil); //left-hand side term, DIMENSION(1:NSOIL)
			FortDoubleArray CI = new FortDoubleArray(1, Driver.NSoil); //left-hand side term, DIMENSION(1:NSOIL)

			double FFF; //runoff decay factor (m-1)
			double RSBMX; //baseflow coefficient [mm/s]

			double PDDUM; //infiltration rate at surface (m/s)
			double FICE; //ice fraction in frozen soil
			double WPLUS; //saturation excess of the total soil [m]
			double RSAT; //accumulation of WPLUS (saturation excess) [m]
			double SICEMAX; //maximum soil ice content (m3/m3)
			double SH2OMIN; //minimum soil liquid water content (m3/m3)
//			double WTSUB; //sum of WCND[K]*DZSNSO[K]
			//double MH2O; //water mass removal (mm)
			double FSAT; //fractional saturated area (-)
			FortDoubleArray MLIQ = new FortDoubleArray(1, Driver.NSoil); //, DIMENSION(1:NSOIL)
			//double XS; //
			//double WATMIN; //
			
			
			FortDoubleArray FCR = new FortDoubleArray(1, Driver.NSoil); //impermeable fraction due to frozen soil, DIMENSION(1:NSOIL)
			int NITER; //iteration times soil moisture (-)
			double SMCTOT; //2-m averaged soil moisture (m3/m3)
			double DZTOT; //2-m soil depth (m)
			

			// ----------------------------------------------------------------------
			RUNSRF = 0.0;
			PDDUM = 0.0;
			RSAT = 0.0;
			QDRAIN = 0;
			RUNSUB = 0;
			double SMCMAX = cell.SMCMAX;
			// for the case when snowmelt water is too large
			//初始化SH2O为低于EPORE的值
//			double EPORE = 0; //effective porosity [m3/m3]
			for (K = 1; K <= Driver.NSoil; K++) {
				double EPORE = Math.Max(1E-4, (SMCMAX - SICE[K]));   //有效空隙率=饱和含水率-土壤冰含量
				RSAT += Math.Max(0, SH2O[K] - EPORE) * DZSNSO[K];  //累积超渗量=(土壤含水率-空隙率)*土层深度，即它是从SH2O中流出的部分
				SH2O[K] = Math.Min(EPORE, SH2O[K]);  //更新SH2O，当原有SH2O超过EPORE时，取EPORE，因为剩余的部分已从SH2O中流出；没超过时按照原有的SH2O，相当于不更新
			}

			//impermeable fraction due to frozen soil
			double A = 4.0;
			double expA = Math.Exp(-A);
			for (K = 1; K <= Driver.NSoil; K++) {
				FICE = Math.Min(1.0, SICE[K] / SMCMAX);  //计算冰含量
				FCR[K] = Math.Max(0.0, Math.Exp(-A * (1 - FICE)) - expA) / (1.0 - expA);  //计算冻土导致的不透水面积占比
			}

			// maximum soil ice content and minimum liquid water of all layers

			SICEMAX = 0.0;
			FCRMAX = 0.0;
			SH2OMIN = SMCMAX;
			for (K = 1; K <= Driver.NSoil; K++) {
				if (SICE[K] > SICEMAX)
					SICEMAX = SICE[K];   //最大冰含量
				if (FCR[K] > FCRMAX)
					FCRMAX = FCR[K];     //最大不透水占比
				if (SH2O[K] < SH2OMIN)
					SH2OMIN = SH2O[K];   //查找到的最小液态水量
			}

			//subsurface runoff for runoff scheme option 2

			if (NoahMP.OPT_RUN == 2) {
				FFF = 2.0;    //径流衰减系数
				RSBMX = 4.0;  //基流系数
				double oldZWT = ZWT;  //计算水位
				ZWT = ZWTEQ(cell, Driver.NSoil, Driver.NSnow, ZSOIL, DZSNSO, SH2O);
				RUNSUB = (1.0 - FCRMAX) * RSBMX * Math.Exp(-NoahMP.TIMEAN) * Math.Exp(-FFF * ZWT);   //地下径流产流 mm/s
				
			}
			//surface runoff and infiltration rate using different schemes

			//jref impermable surface at urban
			if (VEGTYP == Driver.vegparams.ISURBAN)
				FCR[1] = 0.95;

			if (NoahMP.OPT_RUN == 1) {
				FFF = 6.0;  //saturated area decay factor, m-1
				FSAT = NoahMP.FSATMX * Math.Exp(-0.5 * FFF * (ZWT - 2.0));  //计算饱和面积比
				if (QINSUR > 0) {   //当地表输入总水量不为0时
					RUNSRF = QINSUR * ((1.0 - FCR[1]) * FSAT + FCR[1]);
					//PDDUM:infiltration rate at surface (m/s)
					PDDUM = QINSUR - RUNSRF;  //地表向下渗透量PDDUM= m/s
				}
			}

			if (NoahMP.OPT_RUN == 5) {
				FFF = 6.0;
				FSAT = NoahMP.FSATMX * Math.Exp(-0.5 * FFF * Math.Max(-2.0 - ZWT, 0));
				if (QINSUR > 0) {
					RUNSRF = QINSUR * ((1.0 - FCR[1]) * FSAT + FCR[1]);
					PDDUM = QINSUR - RUNSRF;  // m/s
				}
			}
			if (NoahMP.OPT_RUN == 2) {
				FFF = 2.0;
				FSAT = NoahMP.FSATMX * Math.Exp(-0.5 * FFF * ZWT);
				if (QINSUR > 0) {
					RUNSRF = QINSUR * ((1.0 - FCR[1]) * FSAT + FCR[1]);  //前面已算过地下产流量，这里计算地表下渗量
					PDDUM = QINSUR - RUNSRF;                          //地表向下渗量 m/s
				}
			}
			if (NoahMP.OPT_RUN == 3) {
				//不更新SH2O
				INFIL(cell, NSOIL, DT, ZSOIL, SH2O, SICE,
					SICEMAX, QINSUR, out PDDUM, out RUNSRF);
//				if (double.IsNaN(PDDUM))
//					throw new Exception("");
			}


			if (NoahMP.OPT_RUN == 4) {
				SMCTOT = 0;
				DZTOT = 0;
				for (K = 1; K <= Driver.NSoil; K++) {
					DZTOT += DZSNSO[K];
					SMCTOT += SMC[K] * DZSNSO[K];
					if (DZTOT >= 2.0)
						break;
				}
				SMCTOT = SMCTOT / DZTOT;
				FSAT = Math.Pow(Math.Max(0.01, SMCTOT / SMCMAX), 4);        //BATS

				if (QINSUR > 0) {
					RUNSRF = QINSUR * ((1.0 - FCR[1]) * FSAT + FCR[1]);
					PDDUM = QINSUR - RUNSRF;                       // m/s
				}
			}
			// determine iteration times and finer time step

			NITER = 4;
			if (NoahMP.OPT_INF == 1) {     //OPT_INF =2 may cause water imbalance
				NITER = 4;
				if (PDDUM * DT > DZSNSO[1] * SMCMAX)
					NITER *= 2;				
			}

			DTFINE = DT / NITER;			
			double total = TotalWater(PDDUM, DT, ETRANI, QSEVA, SH2O, RSAT, DZSNSO);			
			
			
//			Console.WriteLine("PDDUM=" + PDDUM.ToString("0.00000000000000"));
			//Console.WriteLine("SH2O_sum=" + sum0);
//			Console.WriteLine("QSEVA="+(QSEVA).ToString("0.00000000000000")+" "+"ETRANI="+sum1.ToString("0.00000000000000"));
			//到这里为止，主要误差就在QINSUR上
			double QDRAIN_SAVE = 0.0; //用于累计QDRAIN的储存器
			for (ITER = 0; ITER < NITER; ITER++) {
				//SRT运算已测试通过
				SRT(cell, Driver.NSoil, ZSOIL, DTFINE, PDDUM, ETRANI,
					QSEVA, SH2O, SMC, ZWT, FCR,
					SICEMAX, FCRMAX, ILOC, JLOC, SMCWTD,
					RHSTT, AI, BI, CI, out QDRAIN,
					WCND);

//				if (QDRAIN > 1e-7)
//					throw new Exception();

				//会更新SH2O
				SSTEP(cell, Driver.NSoil, Driver.NSnow, DTFINE, ZSOIL, DZSNSO,
					SICE, ILOC, JLOC, ZWT,
					SH2O, AI, BI, CI,
					RHSTT, ref SMCWTD, ref QDRAIN, ref DEEPRECH,
					out WPLUS);
//				if (QDRAIN > 1e-7)
//					throw new Exception();

				RSAT += WPLUS; //累积饱和多余的水
				QDRAIN_SAVE += QDRAIN;
			}
			QDRAIN = QDRAIN_SAVE / NITER;
			//-------------------------------------------------------
			
			FortDoubleArray etrani = new FortDoubleArray(1, Driver.NSoil);
			double total2 = TotalWater(0, DT, etrani, 0, SH2O, RSAT + QDRAIN, DZSNSO);
			
//			Console.WriteLine("Total0=" + total.ToString("0.00000000000000") + "\nTotal1=" + total2.ToString("0.00000000000000"));
			double error = total2 - total;
//			if (Math.Abs(error) / total > 0.5)
//				throw new Exception();
			//Console.WriteLine("error="+error.ToString("0.00000000000000"));
		
			//下面这个校正是完全不对的，原因未知
//			RSAT -= error;//如果计算出来的水产生了误差造成的多余，则从RSAT中扣除，RSAT中不足时，从土壤水中扣除
//			if (RSAT < -1e-10) {				
//				SH2O[1] += RSAT / DZSNSO[1];
//				RSAT = 0;
//			}
//			double total3= TotalWater(0,DT,etrani,0,SH2O,RSAT+QDRAIN,DZSNSO);
//			error = total3 - total;	
			//-----------------------------------------------------
			RUNSRF = RUNSRF * 1000 + RSAT * 1000 / DT;  // m/s -> mm/s
//			if (RUNSRF > 100)
//				throw new Exception("");
			QDRAIN *= 1000; //QDRAIN最终会被加到RUNSUB中
//			if (QDRAIN > 1e-4)
//					throw new Exception();


			if (NoahMP.OPT_RUN == 2) {
				double WTSUB = 0;//sum of WCND[K]*DZSNSO[K]
				for (K = 1; K <= Driver.NSoil; K++) {
					WTSUB += WCND[K] * DZSNSO[K];   //WCND为水力传导度
				}

				for (K = 1; K <= Driver.NSoil; K++) {
					double MH2O = RUNSUB * DT * (WCND[K] * DZSNSO[K]) / WTSUB;       //water mass removal (mm)
					SH2O[K] -= MH2O / (DZSNSO[K] * 1000);  //土壤水去除渗出量
				}
			}

			// Limit MLIQ to be greater than or equal to watmin.
			// Get water needed to bring MLIQ equal WATMIN from lower layer.

			if (NoahMP.OPT_RUN != 1) {
				for (IZ = 1; IZ <= Driver.NSoil; IZ++) {
					MLIQ[IZ] = SH2O[IZ] * DZSNSO[IZ] * 1000;  //计算各层土壤中的液态水量 mm
					if (double.IsNaN(MLIQ[IZ]))
						throw new Exception("Wrong");
				}

				double XS = 0; //从RUNSUB中进入土壤中的一部分水
				double WATMIN = 0.01;           // mm
				for (IZ = 1; IZ <= Driver.NSoil - 1; IZ++) {
					XS = 0;
					if (MLIQ[IZ] < 0)    //这句意味着SH2O的值可能是负的
						XS = WATMIN - MLIQ[IZ];
					
					MLIQ[IZ] += XS;
					MLIQ[IZ + 1] -= XS;
				}

				IZ = Driver.NSoil;
				XS = 0;
				if (MLIQ[IZ] < WATMIN)
					XS = WATMIN - MLIQ[IZ];
				
				
				MLIQ[IZ] += XS;  //加到土壤水总量中
				RUNSUB -= XS / DT;  //更新RUNSUB
				if (NoahMP.OPT_RUN == 5)
					DEEPRECH -= XS * 1E-3;
				for (IZ = 1; IZ <= Driver.NSoil; IZ++) {
					SH2O[IZ] = MLIQ[IZ] / (DZSNSO[IZ] * 1000);
					if (double.IsNaN(SH2O[IZ]))
						throw new Exception("Wrong");
				}				
			}
			double ttt = RUNSRF * 10800;
		}
				/// <summary>
		/// 计算水力传导度，calculate soil water diffusivity and soil hydraulic conductivity.
		/// </summary>
		/// <param name="WDF"></param>
		/// <param name="WCND"></param>
		/// <param name="SMC"></param>
		/// <param name="SICE"></param>
		public static void WDFCND2(out double WDF, out double WCND, double SMC, double SICE, GridCell cell)
		{
			// soil water diffusivity
			double FACTR = Math.Max(0.01, SMC / cell.SMCMAX);
			double	EXPON = cell.BEXP + 2.0;
			WDF = cell.DWSAT * Math.Pow(FACTR, EXPON);
			double VKWGT = 0;
			if (SICE > 0) {
				VKWGT = 1 / (1 + Math.Pow(500 * SICE, 3));
				WDF = VKWGT * WDF + (1 - VKWGT) * cell.DWSAT * Math.Pow(0.2 / cell.SMCMAX, EXPON);
			}

			// hydraulic conductivity
			EXPON = 2.0 * cell.BEXP + 3.0;
			WCND = cell.DKSAT * Math.Pow(FACTR, EXPON);
		}
		/// <summary>
		///  compute inflitration rate at soil surface and surface runoff
		/// </summary>
		/// <param name="cell"></param>
		/// <param name="NSOIL">no. of soil layers</param>
		/// <param name="DT"></param>
		/// <param name="ZSOIL">depth of soil layer-bottom [m]</param>
		/// <param name="SH2O">soil liquid water content [m3/m3]</param>
		/// <param name="SICE"></param>
		/// <param name="SICEMAX">maximum soil ice content (m3/m3)</param>
		/// <param name="QINSUR">water input on soil surface [mm/s]</param>
		/// <param name="PDDUM">infiltration rate at surface</param>
		/// <param name="RUNSRF">surface runoff [mm/s] </param>
		public static	void INFIL(GridCell cell, int NSOIL, double DT, FortDoubleArray ZSOIL, FortDoubleArray SH2O, FortDoubleArray SICE,  //in
			double SICEMAX, double QINSUR,                          //in
			out double PDDUM, out double RUNSRF)                           //out
		{

// inputs
//  INTEGER,                  INTENT(IN) :: NSOIL  //no. of soil layers
//  REAL,                     INTENT(IN) :: DT     //time step (sec)
//  REAL, DIMENSION(1:NSOIL), INTENT(IN) :: ZSOIL  //depth of soil layer-bottom [m]
//  REAL, DIMENSION(1:NSOIL), INTENT(IN) :: SH2O   //soil liquid water content [m3/m3]
//  REAL, DIMENSION(1:NSOIL), INTENT(IN) :: SICE   //soil ice content [m3/m3]
//  REAL,                     INTENT(IN) :: QINSUR //water input on soil surface [mm/s]
//  REAL,                     INTENT(IN) :: SICEMAX//maximum soil ice content (m3/m3)
//
//// outputs
//  REAL,                    INTENT(OUT) :: RUNSRF //surface runoff [mm/s] 
//  REAL,                    INTENT(OUT) :: PDDUM  //infiltration rate at surface
//
//// locals
//  INTEGER :: IALP1, J, JJ,  K
//  REAL                     :: VAL
//  REAL                     :: DDT
//  REAL                     :: PX
//  REAL                     :: DT1, DD, DICE
//  REAL                     :: FCR
//  REAL                     :: SUM
//  REAL                     :: ACRT
//  REAL                     :: WDF
//  REAL                     :: WCND
//  REAL                     :: SMCAV
//  REAL                     :: INFMAX
			
			int CVFRZ = 3;
			
			if (QINSUR <= 1e-10) {  
				RUNSRF = 0;
				PDDUM = 0;
				return;
			}
			double DT1 = DT / 86400;
			double SMCAV = cell.SMCMAX - cell.SMCWLT;

			// maximum infiltration rate
			FortDoubleArray DMAX = new  FortDoubleArray(1, NSOIL);
			//先计算第1层
			DMAX[1] = -ZSOIL[1] * SMCAV;
			double DICE = -ZSOIL[1] * SICE[1];
			DMAX[1] *= (1.0 - (SH2O[1] + SICE[1] - cell.SMCWLT) / SMCAV);

			double DD = DMAX[1];

			//再计算其它层
			for (int K = 2; K <= NSOIL; K++) {
				DICE += (ZSOIL[K - 1] - ZSOIL[K]) * SICE[K];
				DMAX[K] = (ZSOIL[K - 1] - ZSOIL[K]) * SMCAV;
				DMAX[K] *= (1.0 - (SH2O[K] + SICE[K] - cell.SMCWLT) / SMCAV);
				DD += DMAX[K];
			}
			
			//Arnault et al. (2015) and Yucel et al. (2015):参数REFKDT影响KDT参数。它是影响地表产流的最敏感参数
			double VAL = (1 - Math.Exp(-cell.KDT * DT1));
			double DDT = DD * VAL;
			double PX = Math.Max(0, QINSUR * DT);
			double INFMAX = (PX * (DDT / (PX + DDT))) / DT;

// impermeable fraction due to frozen soil

			double FCR = 1;
			if (DICE > 1E-2) {
				double ACRT = CVFRZ * cell.FRZX / DICE;
				double SUM = 1;
				double	IALP1 = CVFRZ - 1;
				for (int J = 1; J <= IALP1; J++) {
					int	K = 1;
					for (int JJ = J + 1; JJ <= IALP1; JJ++) {
						K *= JJ;
					}
					SUM += Math.Pow(ACRT, CVFRZ - J) / K;
				}
				FCR = 1 - Math.Exp(-ACRT) * SUM;
			}

// correction of infiltration limitation

			INFMAX *= FCR;
// jref for urban areas
//       if (VEGTYP == ISURBAN ) INFMAX == INFMAX * 0.05

			double WDF = 0;
			double WCND = 0;
			WDFCND2(out WDF, out WCND, SH2O[1], SICEMAX, cell);
			INFMAX = Math.Max(INFMAX, WCND);
			INFMAX = Math.Min(INFMAX, PX);

			RUNSRF = Math.Max(0, QINSUR - INFMAX);
			PDDUM = QINSUR - RUNSRF;
		}

		static void CALHUM(double SFCTMP, double SFCPRS, out double Q2SAT, out double DQSDT2)
		{

			//        IMPLICIT NONE
//
			//        REAL, INTENT(IN)       :: SFCTMP, SFCPRS
			//        REAL, INTENT(OUT)      :: Q2SAT, DQSDT2
			double A2 = 17.67;
			double A4 = 29.65;

			double EPSILON = 0.622;
			double E0 = 0.611;
			double A3 = 273.15;
			double A23M4 = A2 * (A3 - A4);
			double RV = 461;
			double ELWV = 2.501E6;
			//        REAL                   :: ES, SFCPRSX

			// Q2SAT: saturated mixing ratio
			double ES = E0 * Math.Exp(ELWV / RV * (1 / A3 - 1 / SFCTMP));
			// convert SFCPRS from Pa to KPa
			double SFCPRSX = SFCPRS * 0.001;
			Q2SAT = EPSILON * ES / (SFCPRSX - ES);
			// convert from  g/g to g/kg
			Q2SAT *= 1000;
			// Q2SAT is currently a 'mixing ratio'

			// DQSDT2 is calculated assuming Q2SAT is a specific humidity
			DQSDT2 = (Q2SAT / (1 + Q2SAT)) * A23M4 / Math.Pow(SFCTMP - A4, 2);

			// DG Q2SAT needs to be in g/g when returned for SFLX
			Q2SAT /= 1000;

		}
		/// <summary>
		/// 经测试，该函数结果已基本正确，仅RHSTT输出有细微差别。
		/// 这里面SH2O会被更新. AI,BI,CI是在SRT中生成，在这里做了修改
		/// </summary>
		/// <param name="NSOIL"></param>
		/// <param name="NSNOW"></param>
		/// <param name="DT"></param>
		/// <param name="ZSOIL"></param>
		/// <param name="DZSNSO"></param>
		/// <param name="SICE"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="ZWT"></param>
		/// <param name="SH2O"></param>
		/// <param name="AI"></param>
		/// <param name="BI"></param>
		/// <param name="CI"></param>
		/// <param name="RHSTT"></param>
		/// <param name="SMCWTD"></param>
		/// <param name="QDRAIN"></param>
		/// <param name="DEEPRECH"></param>
		/// <param name="WPLUS"></param>
		public static	void SSTEP(GridCell cell, int NSOIL, int NSNOW, double DT, FortDoubleArray ZSOIL, FortDoubleArray DZSNSO,
			FortDoubleArray SICE, int ILOC, int JLOC, double ZWT,         //in
			FortDoubleArray SH2O, FortDoubleArray AI, FortDoubleArray BI, FortDoubleArray CI,  //inout
			FortDoubleArray RHSTT, ref double SMCWTD, ref double QDRAIN, ref double DEEPRECH,
			out double WPLUS)
		{
			// ----------------------------------------------------------------------
			// 
			/*/input
				INTEGER,                         INTENT(IN) :: ILOC   //grid index
				INTEGER,                         INTENT(IN) :: JLOC   //grid index
				INTEGER,                         INTENT(IN) :: NSOIL  //
				INTEGER,                         INTENT(IN) :: NSNOW  //
				REAL, INTENT(IN)                            :: DT
				REAL, INTENT(IN)                            :: ZWT
				REAL, DIMENSION(       1:NSOIL), INTENT(IN) :: ZSOIL
				REAL, DIMENSION(       1:NSOIL), INTENT(IN) :: SICE
				REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DZSNSO // snow/soil layer thickness [m]

				//input and output
				REAL, DIMENSION(1:NSOIL), INTENT(INOUT) :: SH2O
				REAL, DIMENSION(1:NSOIL), INTENT(INOUT) :: SMC
				REAL, DIMENSION(1:NSOIL), INTENT(INOUT) :: AI
				REAL, DIMENSION(1:NSOIL), INTENT(INOUT) :: BI
				REAL, DIMENSION(1:NSOIL), INTENT(INOUT) :: CI
				REAL, DIMENSION(1:NSOIL), INTENT(INOUT) :: RHSTT
				REAL                    , INTENT(INOUT) :: SMCWTD
				REAL                    , INTENT(INOUT) :: QDRAIN
				REAL                    , INTENT(INOUT) :: DEEPRECH
				//output
				REAL, INTENT(OUT)                       :: WPLUS     //saturation excess water (m)
				*/
			//local
			int K;
			FortDoubleArray RHSTTIN = new FortDoubleArray(1, Driver.NSoil);
			FortDoubleArray CIIN = new FortDoubleArray(1, Driver.NSoil);
			double SMCMAX = cell.SMCMAX;

			
			double WMINUS;
			// ----------------------------------------------------------------------
			WPLUS = 0.0;
			for (K = 1; K <= NSOIL; K++) {
				RHSTT[K] *= DT;
				AI[K] *= DT;
				BI[K] = 1 + BI[K] * DT;
				CI[K] *= DT;
			}

			// copy values for input variables before calling rosr12

			for (K = 1; K <= NSOIL; K++) {
				RHSTTIN[K] = RHSTT[K];
				CIIN[K] = CI[K];
			}

			//ROSR12 to solve the tri-diagonal matrix

			//不考虑雪层，只求解四层土壤
			ROSR12(CI, AI, BI, CIIN, RHSTTIN, RHSTT, 1, NSOIL, 0);
			for (K = 1; K <= NSOIL; K++) {
				SH2O[K] += CI[K];
			}

			//  excessive water above saturation in a layer is moved to
			//  its unsaturated layer like in a bucket

			//gmmwith opt_run=5 there is soil moisture below nsoil, to the water table
			if (NoahMP.OPT_RUN == 5) {

				//update smcwtd
				if (ZWT < ZSOIL[NSOIL] - DZSNSO[NSOIL]) {
					//accumulate qdrain to update deep water table and soil moisture later
					DEEPRECH += DT * QDRAIN;
				} else {
					SMCWTD += DT * QDRAIN / DZSNSO[NSOIL];
					WPLUS = Math.Max((SMCWTD - SMCMAX), 0.0) * DZSNSO[NSOIL];
					WMINUS = Math.Max((1E-4 - SMCWTD), 0.0) * DZSNSO[NSOIL];

					SMCWTD = Math.Max(Math.Min(SMCWTD, SMCMAX), 1E-4);
					SH2O[NSOIL] += WPLUS / DZSNSO[NSOIL];
					//reduce fluxes at the bottom boundaries accordingly
					QDRAIN -= WPLUS / DT;
					DEEPRECH -= WMINUS;					
				}

			}
			double EPORE = 0;
			for (K = NSOIL; K >= 2; K--) {
				EPORE = Math.Max(1E-4, (SMCMAX - SICE[K]));    //饱和土壤含水率减去土壤冰，得到空隙水容量,SH2O的最终值不能大于EPORE
				WPLUS = Math.Max(SH2O[K] - EPORE, 0.0) * DZSNSO[K];  //土壤水含量减去孔隙水量，就是多余的水，多余出来的水要往上面的土壤层分配
				SH2O[K] = Math.Min(EPORE, SH2O[K]);          //土壤水不能低于空隙水含量
				SH2O[K - 1] += WPLUS / DZSNSO[K - 1];        //多出的水分配至上层土壤中
			}
			//处理最顶层土壤里的水
			EPORE = Math.Max(1E-4, (SMCMAX - SICE[1]));
			WPLUS = Math.Max((SH2O[1] - EPORE), 0.0) * DZSNSO[1];  //最顶层多余的水将作为径流
			SH2O[1] = Math.Min(EPORE, SH2O[1]);			
			
		}
		/// <summary>
		/// 计算气孔阻抗和光合作用量
		/// </summary>
		/// <param name="PAR"></param>
		/// <param name="SFCTMP"></param>
		/// <param name="RCSOIL"></param>
		/// <param name="EAH"></param>
		/// <param name="SFCPRS"></param>
		/// <param name="RC"></param>
		/// <param name="PSN"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		public static void CANRES(GridCell cell, double PAR, double SFCTMP, double RCSOIL, double EAH, double SFCPRS,
			out   double  RC, out double PSN, int ILOC, int JLOC)           //out
		{
			// --------------------------------------------------------------------------------------------------
			// 
			// --------------------------------------------------------------------------------------------------
			// source:  Jarvis (1976), Noilhan and Planton (1989, MWR), Jacquemin and
			// Noilhan (1990, BLM). Chen et al (1996, JGR, Vol 101(D3), 7251-7268),
			// eqns 12-14 and table 2 of sec. 3.1.2
			// --------------------------------------------------------------------------------------------------
			//niu    USE module_Noahlsm_utility
			// --------------------------------------------------------------------------------------------------
			// --------------------------------------------------------------------------------------------------
			// inputs
//
			//    INTEGER,                  INTENT(IN)  :: ILOC   //grid index
			//    INTEGER,                  INTENT(IN)  :: JLOC   //grid index
			//    REAL,                     INTENT(IN)  :: PAR    //par absorbed per unit sunlit lai (w/m2)
			//    REAL,                     INTENT(IN)  :: SFCTMP //canopy air temperature
			//    REAL,                     INTENT(IN)  :: SFCPRS //surface pressure (pa)
			//    REAL,                     INTENT(IN)  :: EAH    //water vapor pressure (pa)
			//    REAL,                     INTENT(IN)  :: RCSOIL //soil moisture stress factor
//
			//outputs
			//  REAL,                     INTENT(OUT) :: RC     //canopy resistance per unit LAI
			//    REAL,                     INTENT(OUT) :: PSN    //foliage photosynthesis (umolco2/m2/s)
//
			//local
//
			//    REAL                                  :: RCQ
			//    REAL                                  :: RCS
			//    REAL                                  :: RCT
			//    REAL                                  :: FF
			//    REAL                                  :: Q2     //water vapor mixing ratio (kg/kg)
			double Q2SAT = 0; //saturation Q2
			double DQSDT2 = 0; //d(Q2SAT)/d(T)

			// RSMIN, RSMAX, TOPT, RGL, HS are canopy stress parameters set in REDPRM
			// ----------------------------------------------------------------------
			// initialize canopy resistance multiplier terms.
			// ----------------------------------------------------------------------
			RC = 0.0;
			double RCS = 0.0;
			double RCT = 0.0;
			double RCQ = 0.0;

			//  compute Q2 and Q2SAT

			double Q2 = 0.622 * EAH / (SFCPRS - 0.378 * EAH); //specific humidity [kg/kg]
			Q2 = Q2 / (1.0 + Q2);                        //mixing ratio [kg/kg]

			CALHUM(SFCTMP, SFCPRS, out Q2SAT, out DQSDT2);

			// contribution due to incoming solar radiation

			double FF = 2.0 * PAR / cell.RGL;
			RCS = (FF + cell.RSMIN / REDPRM.RSMAX) / (1.0 + FF);
			RCS = Math.Max(RCS, 0.0001);

			// contribution due to air temperature

			RCT = 1.0 - 0.0016 * (Math.Pow(REDPRM.TOPT - SFCTMP, 2.0));
			RCT = Math.Max(RCT, 0.0001);

			// contribution due to vapor pressure deficit

			RCQ = 1.0 / (1.0 + cell.HS * Math.Max(0, Q2SAT - Q2));
			RCQ = Math.Max(RCQ, 0.01);

			// determine canopy resistance due to all factors

			RC = cell.RSMIN / (RCS * RCT * RCQ * RCSOIL);
			PSN = -999.99;       // PSN not applied for dynamic carbon
		}
		/// <summary>
		/// 植被通量
		/// </summary>
		/// <param name="NSNOW"></param>
		/// <param name="NSOIL"></param>
		/// <param name="ISNOW"></param>
		/// <param name="VEGTYP"></param>
		/// <param name="VEG"></param>
		/// <param name="DT"></param>
		/// <param name="SAV"></param>
		/// <param name="SAG"></param>
		/// <param name="LWDN"></param>
		/// <param name="UR"></param>
		/// <param name="UU"></param>
		/// <param name="VV"></param>
		/// <param name="SFCTMP"></param>
		/// <param name="THAIR"></param>
		/// <param name="QAIR"></param>
		/// <param name="EAIR"></param>
		/// <param name="RHOAIR"></param>
		/// <param name="SNOWH"></param>
		/// <param name="VAI"></param>
		/// <param name="GAMMAV"></param>
		/// <param name="GAMMAG"></param>
		/// <param name="FWET"></param>
		/// <param name="LAISUN"></param>
		/// <param name="LAISHA"></param>
		/// <param name="CWP"></param>
		/// <param name="DZSNSO"></param>
		/// <param name="HTOP"></param>
		/// <param name="ZLVL"></param>
		/// <param name="ZPD"></param>
		/// <param name="Z0M"></param>
		/// <param name="FVEG"></param>
		/// <param name="Z0MG"></param>
		/// <param name="EMV"></param>
		/// <param name="EMG"></param>
		/// <param name="CANLIQ"></param>
		/// <param name="CANICE"></param>
		/// <param name="STC"></param>
		/// <param name="DF"></param>
		/// <param name="RSSUN">sunlit leaf stomatal resistance (s/m)</param>
		/// <param name="RSSHA">shaded leaf stomatal resistance (s/m)</param>
		/// <param name="RSURF"></param>
		/// <param name="LATHEAV"></param>
		/// <param name="LATHEAG"></param>
		/// <param name="PARSUN">直射光合有效辐射</param>
		/// <param name="PARSHA">阴影光合有效辐射</param>
		/// <param name="IGS"></param>
		/// <param name="FOLN"></param>
		/// <param name="CO2AIR"></param>
		/// <param name="O2AIR"></param>
		/// <param name="BTRAN"></param>
		/// <param name="SFCPRS"></param>
		/// <param name="RHSUR"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="Q2"></param>
		/// <param name="EAH"></param>
		/// <param name="TAH"></param>
		/// <param name="TV"></param>
		/// <param name="TG"></param>
		/// <param name="CM"></param>
		/// <param name="CH"></param>
		/// <param name="DX"></param>
		/// <param name="DZ8W"></param>
		/// <param name="TAUXV"></param>
		/// <param name="TAUYV"></param>
		/// <param name="IRG"></param>
		/// <param name="IRC"></param>
		/// <param name="SHG"></param>
		/// <param name="SHC"></param>
		/// <param name="EVG"></param>
		/// <param name="EVC"></param>
		/// <param name="TR"></param>
		/// <param name="GHV"></param>
		/// <param name="T2MV"></param>
		/// <param name="PSNSUN"></param>
		/// <param name="PSNSHA"></param>
		/// <param name="QC"></param>
		/// <param name="PBLH"></param>
		/// <param name="QSFC"></param>
		/// <param name="PSFC"></param>
		/// <param name="ISURBAN"></param>
		/// <param name="IZ0TLND"></param>
		/// <param name="Q2V"></param>
		/// <param name="CAH2"></param>
		/// <param name="CHLEAF"></param>
		/// <param name="CHUC"></param>
		static void VEGE_FLUX(GridCell cell, int NSNOW, int NSOIL, int ISNOW, int VEGTYP, bool VEG, //in
			double DT, double SAV, double SAG, double LWDN, double UR, //in
			double UU, double VV, double SFCTMP, double THAIR, double QAIR,  //in
			double  EAIR, double RHOAIR, double SNOWH, double VAI, double GAMMAV, double GAMMAG,   //in
			double  FWET, double LAISUN, double LAISHA, double CWP, FortDoubleArray DZSNSO,  //in
			double   HTOP, double ZLVL, double ZPD, double Z0M, double FVEG,   //in
			double Z0MG, double EMV, double EMG, double CANLIQ,            //in
			double CANICE, FortDoubleArray STC, FortDoubleArray DF, out double RSSUN, out double RSSHA,  //in
			double   RSURF, double LATHEAV, double LATHEAG, double PARSUN, double PARSHA, double IGS,  //in
			double   FOLN, double CO2AIR, double O2AIR, double BTRAN, double SFCPRS,  //in
			double   RHSUR, int ILOC, int JLOC, double Q2,   //in
			ref double  EAH, ref double TAH, ref double TV, ref double TG, ref double CM,  //inout
			ref double  CH, double DX, double DZ8W,                    //
			out double  TAUXV, out double TAUYV, out double IRG, out double IRC, out double SHG,  //out
			out double SHC, out double EVG, out double EVC, out double TR, out double GH, //out
			out double  T2MV, out double PSNSUN, out double PSNSHA,                    //out
			double QC, double PBLH, ref double QSFC, double PSFC, int ISURBAN, //in
			int IZ0TLND, out double Q2V, out double CAH2, out double CHLEAF, out double CHUC)
			//inout
		{
//			if (TV < 260) {
//				TV = 260;
//			}
//
			GH = 0;
			EVG = 0;
			
			double CAH = 0;    //sensible heat conductance, canopy air to ZLVL air (m/s)
			double U10V = 0;    //10 m wind speed in eastward dir (m/s)
			double V10V = 0;    //10 m wind speed in eastward dir (m/s)

			double Z0H = 0;          //roughness length, sensible heat (m)
			double Z0HG = 0;         //roughness length, sensible heat (m)
			double RB = 0;           //bulk leaf boundary layer resistance (s/m)

			double RAMG = 0;         //aerodynamic resistance for momentum (s/m)
			double RAHG = 0;         //aerodynamic resistance for sensible heat (s/m)
			double RAWG = 0;         //aerodynamic resistance for water vapor (s/m)

			double CSH = 0;          //coefficients for sh as function of ts
			double CEV = 0;          //coefficients for ev as function of esat[ts]
			double CGH;     //coefficients for st as function of ts

			double ESTV = 0;         //saturation vapor pressure at tv (pa)
			double ESTG = 0;         //saturation vapor pressure at tg (pa)
			double DESTV = 0;        //d(es)/dt at ts (pa/k)
			double DESTG = 0;        //d(es)/dt at tg (pa/k)
			double ESATW = 0;        //es for water
			double ESATI = 0;        //es for ice
			double DSATW = 0;        //d(es)/dt at tg (pa/k) for water
			double DSATI = 0;        //d(es)/dt at tg (pa/k) for ice
//
			double FM = 0;           //momentum stability correction, weighted by prior iters
			double FH = 0;           //sen heat stability correction, weighted by prior iters
			double FHG = 0;          //sen heat stability correction, ground
			double HCAN = 0;         //canopy height (m) [note: hcan >= z0mg]
//
			//  REAL :: A            //temporary calculation
			//  REAL :: B            //temporary calculation
			double CVH = 0;          //sensible heat conductance, leaf surface to canopy air (m/s)

			double MOZ = 0;          //Monin-Obukhov stability parameter
			double MOZG = 0;         //Monin-Obukhov stability parameter
			//double MOZOLD = 0;       //Monin-Obukhov stability parameter from prior iteration
			double FM2 = 0;          //Monin-Obukhov momentum adjustment at 2m
			double FH2 = 0;          //Monin-Obukhov heat adjustment at 2m
			double CH2 = 0;          //Surface exchange at 2m
			//  REAL :: THSTAR          //Surface exchange at 2m
			//     	 REAL :: THVAIR
			//  REAL :: THAH
			double RAHC2 = 0;        //aerodynamic resistance for sensible heat (s/m)
			double RAWC2 = 0;        //aerodynamic resistance for water vapor (s/m)
			//  REAL, INTENT(OUT):: CAH2         //sensible heat conductance for diagnostics
			double CH2V = 0;         //exchange coefficient for 2m over vegetation.
			double CQ2V = 0;         //exchange coefficient for 2m over vegetation.

			//jref - NITERC test from 5 to 20
			int NITERC = 20;   //number of iterations for surface temperature
			//jref - NITERG test from 3-5
			int NITERG = 5;   //number of iterations for ground temperature

			// Energy released or consumed by snow  frozen soil

			//double TAH = this.TAH;
			// ---------------------------------------------------------------------------------------------

//			IRG = 0;
			double MPE = 1E-6;
			int LITER = 0;
			double FV = 0.1;

			// ---------------------------------------------------------------------------------------------
			// initialization variables that do not depend on stability iteration
			// ---------------------------------------------------------------------------------------------
			double DTV = 0;
			double DTG = 0;
			double MOZSGN = 0;
			double MOZOLD = 0;
			//感热通量
			double HG = 0;
			double H = 0;
			double QFX = 0;
			
			RSSUN = -9999;//自己加的，仅仅是为了语法
			RSSHA = -9999;
			SHG = -9999;
			SHC = -9999;
			EVC = -9999;
			TR = -9999; //transpiration heat flux (w/m2)[+= to atm]
			IRG = -9999;
			IRC = -9999;
			T2MV = -9999;
			PSNSUN = -9999; //sunlit photosynthesis (umolco2/m2/s)
			PSNSHA = -9999; //shaded photosynthesis (umolco2/m2/s)
			Q2V = -9999;
			CAH2 = -9999;


			// convert grid-cell LAI to the fractional vegetated area (FVEG)
			double VAIE = Math.Min(6, VAI / FVEG);
			double LAISUNE = Math.Min(6, LAISUN / FVEG);
			double LAISHAE = Math.Min(6, LAISHA / FVEG);

			// saturation vapor pressure at ground temperature

			double T = TDC(TG);
			
			ESAT(T, out ESATW, out ESATI, out DSATW, out DSATI);
			if (T > 0) {
				ESTG = ESATW;
			} else {
				ESTG = ESATI;
			}

			//jref - consistent surface specific humidity for sfcdif3 and sfcdif4

			QSFC = 0.622 * EAIR / (PSFC - 0.378 * EAIR);

			// canopy height

			HCAN = HTOP;
			double UC = UR * Math.Log(HCAN / Z0M) / Math.Log(ZLVL / Z0M);
			if ((HCAN - ZPD) <= 0) {
				throw new Exception("植被能量错误");
				//          WRITE(message,*) "CRITICAL PROBLEM: HCAN <= ZPD"
				//          wrf_message ( message )
				//          WRITE(message,*) 'i,j point=',ILOC, JLOC
				//          wrf_message ( message )
				//          WRITE(message,*) 'HCAN  =',HCAN
				//          wrf_message ( message )
				//          WRITE(message,*) 'ZPD   =',ZPD
				//          wrf_message ( message )
				//          write (message, *) 'SNOWH =',SNOWH
				// 	          wrf_message ( message )
				//          wrf_error_fatal ( "CRITICAL PROBLEM IN MODULE_SF_NoahMP3LSM:VEGEFLUX" )
			}

			// prepare for longwave rad.

			double AIR = -EMV * (1 + (1 - EMV) * (1 - EMG)) * LWDN - EMV * EMG * NoahMP.SB * Math.Pow(TG, 4);
			double CIR = (2 - EMV * (1 - EMG)) * EMV * NoahMP.SB;
			
			
			// ---------------------------------------------------------------------------------------------
			//  begin stability iteration
			//Console.WriteLine("NITERC:"+NITERC);
			int ITER = 0;
			for (ITER = 1; ITER <= NITERC; ITER++) {
				if (ITER == 1) {
					Z0H = Z0M;
					Z0HG = Z0MG;
				} else {
					Z0H = Z0M;    //* Math.Exp(-CZIL*0.4*258.2*SQRT(FV*Z0M))
					Z0HG = Z0MG;   //* Math.Exp(-CZIL*0.4*258.2*SQRT(FV*Z0MG))
				}

				// aerodyn resistances between heights zlvl and d+z0v

				if (NoahMP.OPT_SFC == 1) {
					//已调试验证，无误，这是原来NoahMP3中的默认选项
					SFCDIF1.SFCDIF(ITER, SFCTMP, RHOAIR, H, QAIR,  //in
						ZLVL, ZPD, Z0M, Z0H, UR,  //in
						MPE, ILOC, JLOC,                  //in
						ref MOZ, ref MOZSGN, ref FM, ref FH, ref FM2, ref FH2,  //inout
						out CM, out CH, ref FV, out CH2);          //out
				}

				double WSTAR = 0;// 0.949328840;
				if (NoahMP.OPT_SFC == 2) {
					//发现原来FV参数标为了out属性，应该是inout属性才是对的，经验证基本输出正确
					SFCDIF2.SFCDIF(ITER, Z0M, TAH, THAIR, UR,   //in
						NoahMP.CZIL, ZLVL, ILOC, JLOC,          //in
						ref CM, ref CH, ref MOZ, ref WSTAR,     //in
						ref FV);                                //out
					// Undo the multiplication by windspeed that SFCDIF2
					// applies to exchange coefficients CH and CM:
					CH = CH / UR;
					CM = CM / UR;
				}

				if (NoahMP.OPT_SFC == 3) {
					//此选项最好不要用，目前还不收敛
					throw new Exception("此选项最好不要用，目前还不收敛");
					SFCDIF3.SFCDIF(ILOC, JLOC, TAH, QSFC, PSFC,  //in
						PBLH, Z0M, Z0MG, VEGTYP, ISURBAN,  //in
						IZ0TLND, UC, ITER, NITERC, SFCTMP,  //in
						THAIR, QAIR, QC, ZLVL,   //in
						SFCPRS, ref FV, ref CM, ref CH, out CH2V,  //inout
						out CQ2V, out MOZ);                             //out
					// Undo the multiplication by windspeed that SFCDIF3
					// applies to exchange coefficients CH and CM:
					CH = CH / UR;
					CM = CM / UR;
					CH2V = CH2V / UR;
				}

				if (NoahMP.OPT_SFC == 4) {
//					if (double.IsNaN(QSFC))
//						throw new Exception("");
//					Console.WriteLine("Iter="+ITER+ "QSFC:"+QSFC);
					//此选项实际上大概不好用
					SFCDIF4.SFCDIF(ILOC, JLOC, UU, VV, SFCTMP,  //in
						SFCPRS, PSFC, PBLH, DX, Z0M,
						TAH, QAIR, ZLVL, IZ0TLND, ref QSFC,
						H, QFX, ref CM, ref CH, ref CH2V,
						ref CQ2V, ref MOZ, ref FV, out U10V, out V10V);
					
					// Undo the multiplication by windspeed that SFCDIF4
					// applies to exchange coefficients CH and CM:
					CH = CH / UR;
					CM = CM / UR;
					CH2V = CH2V / UR;
				}
				double RAMC = Math.Max(1.0, 1.0 / (CM * UR));
				double RAHC = Math.Max(1.0, 1.0 / (CH * UR));
				double RAWC = RAHC;

				if (NoahMP.OPT_SFC == 3 || NoahMP.OPT_SFC == 4) {
					RAHC2 = Math.Max(1.0, 1.0 / (CH2V * UR));
					RAWC2 = RAHC2;
					CAH2 = 1.0 / RAHC2;
					CQ2V = CAH2;
				}
				// aerodyn resistance between heights z0g and d+z0v, RAG, and leaf
				// boundary layer resistance, RB

				
//				if (double.IsInfinity(TV))
//					throw new Exception("");
				RAGRB(ITER, VAIE, RHOAIR, HG, TAH,  //in
					ZPD, Z0MG, Z0HG, HCAN, UC,  //in
					Z0H, FV, CWP, VEGTYP, MPE, TV,  //in
					ref MOZG, ref FHG, ILOC, JLOC,  //inout
					out RAMG, out RAHG, out RAWG, out RB);           //out

//				if (TV>320)
//					throw new Exception("");
				// es and d(es)/dt evaluated at tv

				T = TDC(TV);
				ESAT(T, out ESATW, out ESATI, out DSATW, out DSATI);
				if (T > 0) {
					ESTV = ESATW;
					DESTV = DSATW;
				} else {
					ESTV = ESATI;
					DESTV = DSATI;
				}
				// stomatal resistance
				
				if (ITER == 1) {
					if (NoahMP.OPT_CRS == 1) {        	// Ball-Berry
						//计算阳光直射下的气孔阻抗(即beta)和光合作用量
						STOMATA(VEGTYP, MPE, PARSUN, FOLN, ILOC, JLOC,  //in
							TV, ESTV, EAH, SFCTMP, SFCPRS,  //in
							O2AIR, CO2AIR, IGS, BTRAN, RB,  //in
							out RSSUN, out PSNSUN);                        //out
//						if (double.IsNaN(EAH))
//							throw new Exception();

						//计算阴影散射下的气孔阻抗和光合作用量
						STOMATA(VEGTYP, MPE, PARSHA, FOLN, ILOC, JLOC,  //in
							TV, ESTV, EAH, SFCTMP, SFCPRS,  //in
							O2AIR, CO2AIR, IGS, BTRAN, RB,  //in
							out	RSSHA, out PSNSHA);                      //out
//						Console.WriteLine("PSNSHA="+PSNSHA+" PSNSUN="+PSNSUN);
					}

					if (NoahMP.OPT_CRS == 2) {        	// Jarvis
						if (NoahMP.DVEG == 2) {
							throw new Exception("CANOPY_STOMATAL_RESISTANCE_OPTION must be 1 when DYNAMIC_VEG_OPTION == 2");
						}
						//计算阳光直射下的气孔阻抗和光合作用量
						CANRES(cell, PARSUN, TV, BTRAN, EAH, SFCPRS,  //in
							out RSSUN, out PSNSUN, ILOC, JLOC);          //out

						//计算阴影散射下的气孔阻抗和光合作用量
						CANRES(cell, PARSHA, TV, BTRAN, EAH, SFCPRS,  //in
							out RSSHA, out PSNSHA, ILOC, JLOC);          //out
						
					}
				} //end if

				// prepare for sensible heat flux above veg.

				CAH = 1 / RAHC;
				CVH = 2 * VAIE / RB;
				CGH = 1 / RAHG;
				double COND = CAH + CVH + CGH;
				double ATA = (SFCTMP * CAH + TG * CGH) / COND;
				double BTA = CVH / COND;
				CSH = (1 - BTA) * RHOAIR * NoahMP.CPAIR * CVH;
				// prepare for latent heat flux above veg.

				double CAW = 1.0 / RAWC;
				double	CEW = FWET * VAIE / RB;
				double CTW = (1 - FWET) * (LAISUNE / (RB + RSSUN) + LAISHAE / (RB + RSSHA));
				double CGW = 1.0 / (RAWG + RSURF);
				COND = CAW + CEW + CTW + CGW;
				double AEA = (EAIR * CAW + ESTG * CGW) / COND;
				double BEA = (CEW + CTW) / COND;
				CEV = (1 - BEA) * CEW * RHOAIR * NoahMP.CPAIR / GAMMAV;   // Barlage: change to vegetation v3.6
				double CTR = (1 - BEA) * CTW * RHOAIR * NoahMP.CPAIR / GAMMAV;

				// evaluate surface fluxes with current temperature and solve for dts

				TAH = ATA + BTA * TV;              // canopy air T.
				EAH = AEA + BEA * ESTV;             // canopy air e
				
				IRC = FVEG * (AIR + CIR * Math.Pow(TV, 4));
				SHC = FVEG * RHOAIR * NoahMP.CPAIR * CVH * (TV - TAH);
				EVC = FVEG * RHOAIR * NoahMP.CPAIR * CEW * (ESTV - EAH) / GAMMAV; // Barlage: change to v in v3.6
				
				TR = FVEG * RHOAIR * NoahMP.CPAIR * CTW * (ESTV - EAH) / GAMMAV;
				if (TV > NoahMP.TFRZ) {
					EVC = Math.Min(CANLIQ * LATHEAV / DT, EVC);    // Barlage: add if block for canice in v3.6
				} else {
					EVC = Math.Min(CANICE * LATHEAV / DT, EVC);
				}

				double TV3 = TV * TV * TV;
				double	B = SAV - IRC - SHC - EVC - TR;                         //additional w/m2
				double	A = FVEG * (4.0 * CIR * TV3 + CSH + (CEV + CTR) * DESTV);//volumetric heat capacity
				//关键是看IRC, SHC, CSH
				
				DTV = B / A;

				IRC += FVEG * 4.0 * CIR * TV3 * DTV;
				SHC += FVEG * CSH * DTV;				
				EVC += FVEG * CEV * DESTV * DTV;
//				Console.WriteLine("EVC="+EVC+"IRC="+IRC.ToString("0.000")+" DTV="+DTV);
				TR += FVEG * CTR * DESTV * DTV;
				TV += DTV;
				//        TAH = ATA + BTA*TV               // canopy air T; update here for consistency

				// for computing M-O length in the next iteration
				H = RHOAIR * NoahMP.CPAIR * (TAH - SFCTMP) / RAHC;
				HG = RHOAIR * NoahMP.CPAIR * (TG - TAH) / RAHG;

				// consistent specific humidity from canopy air vapor pressure
				QSFC = (0.622 * EAH) / (SFCPRS - 0.378 * EAH);

				// added moisture flux for sfcdif4
				if (NoahMP.OPT_SFC == 4) {
					QFX = (QSFC - QAIR) * RHOAIR * CAW; //*CPAIR/GAMMAV
				}				
				
//				if (ITER >= 25 && Math.Abs(DTV) <= 0.001 && Math.Abs(LITER) < 1e-10) {
//					break;
//				}
				if (LITER == 1)
					break;
				if (ITER >= 5 && Math.Abs(DTV) <= 0.01 && LITER == 0) {
					LITER = 1;
				}
			} //END DO loop1 // end stability iteration
			
			
			double TV4 = TV * TV * TV * TV;
			AIR = -EMG * (1 - EMV) * LWDN - EMG * EMV * NoahMP.SB * TV4;
			CIR = EMG * NoahMP.SB;
			CSH = RHOAIR * NoahMP.CPAIR / RAHG;
			CEV = RHOAIR * NoahMP.CPAIR / (GAMMAG * (RAWG + RSURF)); // Barlage: change to ground v3.6
			CGH = 2.0 * DF[ISNOW + 1] / DZSNSO[ISNOW + 1]; //

			
			for (ITER = 1; ITER <= NITERG; ITER++) {
				T = TDC(TG);
				ESAT(T, out ESATW, out ESATI, out DSATW, out DSATI);
				if (T > 0) {
					ESTG = ESATW;
					DESTG = DSATW;
				} else {
					ESTG = ESATI;
					DESTG = DSATI;
				}

				IRG = CIR * Math.Pow(TG, 4) + AIR;
				SHG = CSH * (TG - TAH);
				EVG = CEV * (ESTG * RHSUR - EAH);
				GH = CGH * (TG - STC[ISNOW + 1]);
				

				double TG3 = TG * TG * TG;
				double B = SAG - IRG - SHG - EVG - GH;
				double A = 4.0 * CIR * TG3 + CSH + CEV * DESTG + CGH;
				DTG = B / A;

				IRG += 4.0 * CIR * TG3 * DTG;
				SHG += CSH * DTG;
				EVG += CEV * DESTG * DTG;
				GH += CGH * DTG;
				TG += DTG;
			} // end the second loop

			//     TAH = (CAH*SFCTMP + CVH*TV + CGH*TG)/(CAH + CVH + CGH)

			// if snow on ground and TG > TFRZ: reset TG = TFRZ. reevaluate ground fluxes.

			if (NoahMP.OPT_STC == 1) {
				if (SNOWH > 0.05 && TG > NoahMP.TFRZ) {
					TG = NoahMP.TFRZ;
					IRG = CIR * Math.Pow(TG, 4) - EMG * (1 - EMV) * LWDN - EMG * EMV * NoahMP.SB * Math.Pow(TV, 4);
					SHG = CSH * (TG - TAH);
					EVG = CEV * (ESTG * RHSUR - EAH);
					GH = SAG - (IRG + SHG + EVG);
				}
			}

			// wind stresses

			TAUXV = -RHOAIR * CM * UR * UU;
			TAUYV = -RHOAIR * CM * UR * VV;

			// consistent vegetation air temperature and vapor pressure since TG is not consistent with the TAH/EAH
			// calculation.
			//     TAH = SFCTMP + (SHG+SHC)/(RHOAIR*CPAIR*CAH)
			//     TAH = SFCTMP + (SHG*FVEG+SHC)/(RHOAIR*CPAIR*CAH) // ground flux need fveg
			//     EAH = EAIR + (EVC+FVEG*(TR+EVG))/(RHOAIR*CAW*CPAIR/GAMMAG )
			//     QFX = (QSFC-QAIR)*RHOAIR*CAW //*CPAIR/GAMMAG

			// 2m temperature over vegetation ( corrected for low CQ2V values )
			if (NoahMP.OPT_SFC == 1 || NoahMP.OPT_SFC == 2) {
				//      CAH2 = FV*1.0/VKC*LOG((2.0+Z0H)/Z0H)
				CAH2 = FV * NoahMP.VKC / Math.Log((2.0 + Z0H) / Z0H);
				CAH2 = FV * NoahMP.VKC / (Math.Log((2.0 + Z0H) / Z0H) - FH2);
				CQ2V = CAH2;
				if (CAH2 < 1E-5) {
					T2MV = TAH;
					//         Q2V  = (EAH*0.622/(SFCPRS - 0.378*EAH))
					Q2V = QSFC;
				} else {
					T2MV = TAH - (SHG + SHC / FVEG) / (RHOAIR * NoahMP.CPAIR) * 1 / CAH2;
					//         Q2V = (EAH*0.622/(SFCPRS - 0.378*EAH))- QFX/(RHOAIR*FV)* 1.0/VKC * LOG((2.0+Z0H)/Z0H)
					Q2V = QSFC - ((EVC + TR) / FVEG + EVG) / (LATHEAV * RHOAIR) * 1 / CQ2V;
				}
			}

			// myj/ysu consistent 2m temperature over vegetation (if CQ2V < 1e-5? )
			if (NoahMP.OPT_SFC == 3 || NoahMP.OPT_SFC == 4) {
				if (CAH2 < 1E-5) {
					T2MV = TAH;
					Q2V = (EAH * 0.622 / (SFCPRS - 0.378 * EAH));
				} else {
					T2MV = TAH - (SHG + SHC) / (RHOAIR * NoahMP.CPAIR * CAH2);
					Q2V = (EAH * 0.622 / (SFCPRS - 0.378 * EAH)) - QFX / (RHOAIR * CQ2V);
				}
			}

			// update CH for output
			CH = CAH;
			CHLEAF = CVH;
			CHUC = 1 / RAHG;

		}
		
		
		/// <summary>
		/// 已调试过，应该没问题了
		/// </summary>
		/// <param name="NSNOW"></param>
		/// <param name="NSOIL"></param>
		/// <param name="ISNOW"></param>
		/// <param name="DT"></param>
		/// <param name="AI"></param>
		/// <param name="BI"></param>
		/// <param name="CI"></param>
		/// <param name="RHSTS"></param>
		/// <param name="STC"></param>
		public static void HSTEP(int NSNOW, int NSOIL, int ISNOW, double DT,
			FortDoubleArray  AI, FortDoubleArray BI, FortDoubleArray CI, FortDoubleArray RHSTS,
			FortDoubleArray  STC)
		{
			// ----------------------------------------------------------------------
			// CALCULATE/UPDATE THE SOIL TEMPERATURE FIELD.
			// ----------------------------------------------------------------------
			// ----------------------------------------------------------------------
			// input
//
			//    INTEGER,                         INTENT(IN)    :: NSOIL
			//    INTEGER,                         INTENT(IN)    :: NSNOW
			//    INTEGER,                         INTENT(IN)    :: ISNOW
			//    REAL,                            INTENT(IN)    :: DT

			// output  input
			//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: RHSTS
			//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: AI
			//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: BI
			//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: CI
			//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC

			// local
			
			FortDoubleArray RHSTSIN = new FortDoubleArray(1 - NSNOW, NSOIL);
			FortDoubleArray CIIN = new FortDoubleArray(1 - NSNOW, NSOIL);
			
			for (int K = ISNOW + 1; K <= NSOIL; K++) {
				
				RHSTS[K] *= DT;
				AI[K] *= DT;
				BI[K] = 1 + BI[K] * DT;
				CI[K] *= DT;
			}
			// copy values for input variables beforeto rosr12
			for (int K = ISNOW + 1; K <= NSOIL; K++) {
				
				RHSTSIN[K] = RHSTS[K];
				CIIN[K] = CI[K];
			}
			// solve the tri-diagonal matrix equation
			//考虑了雪层，求解7层雪+土壤，但是从顶层的ISNOW起步，至最底层结束
			//ROSR12(CI, AI, BI, CIIN, RHSTSIN, RHSTS, ISNOW + 1, NSNOW);
			ROSR12(CI, AI, BI, CIIN, RHSTSIN, RHSTS, ISNOW + 1, NSOIL, NSNOW);
			// update snow & soil temperature
			for (int K = ISNOW + 1; K <= NSOIL; K++) {
				STC[K] += CI[K];
			}
		}
		/// <summary>
		/// NTOP雪顶或土顶层数，包含了这个数,2022年11月16日检查了一次
		/// </summary>
		/// <param name="P"></param>
		/// <param name="A"></param>
		/// <param name="B"></param>
		/// <param name="C"></param>
		/// <param name="D"></param>
		/// <param name="DELTA"></param>
		/// <param name="NTOP"></param>
		/// <param name="NBottom">按照原有Fortran方式，包含这个数，即NSnow+NSoil-1</param>
		/// <param name="NSNOW"></param>
		static void ROSR12(FortDoubleArray P, FortDoubleArray A, FortDoubleArray B, FortDoubleArray C, FortDoubleArray D, FortDoubleArray DELTA, int NTOP, int NSOIL, int NSNOW)
		{
			// ----------------------------------------------------------------------
			//void ROSR12
			// ----------------------------------------------------------------------
			// INVERT (SOLVE) THE TRI-DIAGONAL MATRIX PROBLEM SHOWN BELOW:
			// ###                                            ### ###  ###   ###  ###
			// #B[1], C[1],  0  ,  0  ,  0  ,   . . .  ,    0   # #      #   #      #
			// #A(2), B(2), C(2),  0  ,  0  ,   . . .  ,    0   # #      #   #      #
			// # 0  , A(3), B(3), C(3),  0  ,   . . .  ,    0   # #      #   # D(3) #
			// # 0  ,  0  , A(4), B(4), C(4),   . . .  ,    0   # # P(4) #   # D(4) #
			// # 0  ,  0  ,  0  , A(5), B(5),   . . .  ,    0   # # P(5) #   # D(5) #
			// # .                                          .   # #  .   # = #   .  #
			// # .                                          .   # #  .   #   #   .  #
			// # .                                          .   # #  .   #   #   .  #
			// # 0  , . . . , 0 , A(M-2), B(M-2), C(M-2),   0   # #P(M-2)#   #D(M-2)#
			// # 0  , . . . , 0 ,   0   , A(M-1), B(M-1), C(M-1)# #P(M-1)#   #D(M-1)#
			// # 0  , . . . , 0 ,   0   ,   0   ,  A(M) ,  B(M) # # P(M) #   # D(M) #
			// ###                                            ### ###  ###   ###  ###
			// ----------------------------------------------------------------------
			

			//    INTEGER, INTENT(IN)   :: NTOP
			//    INTEGER, INTENT(IN)   :: NSOIL,NSNOW
			//    INTEGER               :: K, KK
			//    REAL, DIMENSION(-NSNOW+1:NSOIL),INTENT(IN):: A, B, D
			//    REAL, DIMENSION(-NSNOW+1:NSOIL),INTENT(INOUT):: C,P,DELTA
			
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

			for (int K = NTOP + 1; K <= NSOIL; K++) {
				double v = 1.0 / (B[K] + A[K] * P[K - 1]);
//				if(v==0)
//					throw new Exception("");
				P[K] = -C[K] * v;
				DELTA[K] = (D[K] - A[K] * DELTA[K - 1]) * v;
			}
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
			
		}
		/// <summary>
		/// 雪和土壤温度
		/// </summary>
		/// <param name="ICE"></param>
		/// <param name="NSOIL"></param>
		/// <param name="NSNOW"></param>
		/// <param name="ISNOW"></param>
		/// <param name="IST"></param>
		/// <param name="TBOT"></param>
		/// <param name="ZSNSO"></param>
		/// <param name="SSOIL"></param>
		/// <param name="DF"></param>
		/// <param name="HCPCT"></param>
		/// <param name="ZBOT"></param>
		/// <param name="SAG"></param>
		/// <param name="DT"></param>
		/// <param name="SNOWH"></param>
		/// <param name="DZSNSO"></param>
		/// <param name="TG"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="STC"></param>
		static void TSNOSOI(int ICE, int NSOIL, int NSNOW, int ISNOW, double IST, double  //in
		             TBOT, FortDoubleArray ZSNSO, double SSOIL, FortDoubleArray DF, FortDoubleArray HCPCT, double  //in
		             ZBOT, double SAG, double DT, double SNOWH, FortDoubleArray DZSNSO, double  //in
		             TG, double ILOC, double JLOC,                    //in
			FortDoubleArray STC)                                       //inout
		{
			//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC  //唯一的输出变量

			int n = Driver.NSnow + Driver.NSoil;

			double EFLXB = 0; //energy influx from soil bottom (w/m2)
			//double PHI=0;   //light through water (w/m2)
//
			FortDoubleArray TBEG = new FortDoubleArray(1 - NSNOW, NSOIL);
			double ERR_EST = 0; //heat storage error  (w/m2)
			double SSOIL2 = 0;  //ground heat flux (w/m2) (for energy check)
			double EFLXB2 = 0;  //heat flux from the bottom (w/m2) (for energy check)
			string message = "";
			// ----------------------------------------------------------------------
			// compute solar penetration through water, needs more work

			//原NoahMP3中此变量不起任何传递参数的作用
			FortDoubleArray PHI = new FortDoubleArray(1 - NSNOW, NSOIL);

			// adjust ZBOT from soil surface to ZBOTSNO from snow surface

			double ZBOTSNO = ZBOT - SNOWH;    //from snow surface

			// snow/soil heat storage for energy balance check
			for (int IZ = ISNOW + 1; IZ <= Driver.NSoil; IZ++) {
				TBEG[IZ] = STC[IZ];
			}
			// compute soil temperatures

			FortDoubleArray AI = new FortDoubleArray(1 - NSNOW, NSOIL);
			FortDoubleArray BI = new FortDoubleArray(1 - NSNOW, NSOIL);
			FortDoubleArray CI = new FortDoubleArray(1 - NSNOW, NSOIL);
			FortDoubleArray RHSTS = new FortDoubleArray(1 - NSNOW, NSOIL);
			HRT(Driver.NSnow, Driver.NSoil, ISNOW, ZSNSO,
				STC, TBOT, ZBOTSNO, DT,
				DF, HCPCT, SSOIL, PHI,
				AI, BI, CI, RHSTS,
				out EFLXB);
			//if (double.IsNaN(STC[0]))
//					throw new Exception("");
			HSTEP(Driver.NSnow, Driver.NSoil, ISNOW, DT, AI, BI, CI, RHSTS, STC);
			//if (double.IsNaN(STC[0]))
//	throw new Exception("");
			// update ground heat flux just for energy check, but not for final output
			// otherwise, it would break the surface energy balance

			if (NoahMP.OPT_TBOT == 1) {
				EFLXB2 = 0;
			} else if (NoahMP.OPT_TBOT == 2) {
				EFLXB2 = DF[Driver.NSoil] * (TBOT - STC[Driver.NSoil]) / (0.5 * (ZSNSO[Driver.NSoil - 1] + ZSNSO[Driver.NSoil]) - ZBOTSNO);
			}
			
		}
		/// <summary>
		/// 有关SNOW WATER这部分，因为在测试时，只考虑了ISNOW=3的情况， 即实际上不包含雪的问题。所以需要今后去针对ISNOW<3的情况进行调试
		/// 2022年11月16日又检查了一遍
		/// </summary>
		/// <param name="NSNOW"></param>
		/// <param name="NSOIL"></param>
		/// <param name="IMELT"></param>
		/// <param name="DT"></param>
		/// <param name="ZSOIL"></param>
		/// <param name="SFCTMP"></param>
		/// <param name="SNOWHIN"></param>
		/// <param name="QSNOW"></param>
		/// <param name="QSNFRO"></param>
		/// <param name="QSNSUB"></param>
		/// <param name="QRAIN"></param>
		/// <param name="FICEOLD"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="ISNOW"></param>
		/// <param name="SNOWH"></param>
		/// <param name="SNEQV"></param>
		/// <param name="SNICE"></param>
		/// <param name="SNLIQ"></param>
		/// <param name="SH2O"></param>
		/// <param name="SICE"></param>
		/// <param name="STC"></param>
		/// <param name="ZSNSO"></param>
		/// <param name="DZSNSO"></param>
		/// <param name="QSNBOT"></param>
		/// <param name="SNOFLOW"></param>
		/// <param name="PONDING1"></param>
		/// <param name="PONDING2"></param>
		static void SNOWWATER(int NSNOW, int NSOIL, FortIntArray IMELT, double DT, FortDoubleArray ZSOIL,
			double SFCTMP, double SNOWHIN, double QSNOW, double QSNFRO, double QSNSUB,
			double QRAIN, FortDoubleArray FICEOLD, int ILOC, int JLOC,
			ref int ISNOW, ref double SNOWH, ref double SNEQV, FortDoubleArray SNICE, FortDoubleArray SNLIQ,
			FortDoubleArray SH2O, FortDoubleArray SICE, FortDoubleArray STC, FortDoubleArray ZSNSO, FortDoubleArray DZSNSO,
			out double QSNBOT, out double SNOFLOW, out double PONDING1, out double PONDING2)
		{

			// input
			/*
  int                              ILOC   ;   //grid index
  int                              JLOC   ;   //grid index
  int                              NSNOW  ;   //maximum no. of snow layers
  int                              NSOIL  ;   //no. of soil layers
  intDIMENSION(-NSNOW+1:0) ,      IMELT  ;   //melting state index [0-no melt;1-melt]
  double                                 DT     ;   //time step (s)
  double DIMENSION(       1:NSOIL),      ZSOIL  ;   //depth of layer-bottom from soil surface
  double                                 SFCTMP ;   //surface air temperature [k]
  double                                 SNOWHIN;   //snow depth increasing rate (m/s)
  double                                 QSNOW  ;   //snow at ground srf (mm/s) [+]
  double                                 QSNFRO ;   //snow surface frost rate[mm/s]
  double                                 QSNSUB ;   //snow surface sublimation rate[mm/s]
  double                                 QRAIN  ;   //snow surface rain rate[mm/s]
			double DIMENSION(-NSNOW+1:0)    ,      FICEOLD;   //ice fraction at last timestep

! input & output
  INTEGER,                         INTENT(INOUT) :: ISNOW  !actual no. of snow layers
  REAL,                            INTENT(INOUT) :: SNOWH  !snow height [m]
  REAL,                            INTENT(INOUT) :: SNEQV  !snow water eqv. [mm]
  REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE  !snow layer ice [mm]
  REAL, DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ  !snow layer liquid water [mm]
  REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O   !soil liquid moisture (m3/m3)
  REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SICE   !soil ice moisture (m3/m3)
  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC    !snow layer temperature [k]
  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: ZSNSO  !depth of snow/soil layer-bottom
  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO !snow/soil layer thickness [m]

! output
  REAL,                              INTENT(OUT) :: QSNBOT !melting water out of snow bottom [mm/s]
  REAL,                              INTENT(OUT) :: SNOFLOW!glacier flow [mm]
  REAL,                              INTENT(OUT) :: PONDING1
  REAL,                              INTENT(OUT) :: PONDING2
			 */
			// local
			int IZ, i;
			double BDSNOW;   //bulk density of snow (kg/m3)
			SNOFLOW = 0.0;
			PONDING1 = 0.0;
			PONDING2 = 0.0;
//			if (double.IsNaN(DZSNSO[0]))
//				throw new Exception("");
			SNOWFALL(Driver.NSoil, Driver.NSnow, DT, QSNOW, SNOWHIN,
				SFCTMP, ILOC, JLOC,
				ref ISNOW, ref SNOWH, DZSNSO, STC, SNICE,
				SNLIQ, ref SNEQV);

			if (ISNOW < 0) {
				COMPACT(Driver.NSnow, Driver.NSoil, DT, STC, SNICE,
					SNLIQ, ZSOIL, IMELT, FICEOLD, ILOC, JLOC,
					ISNOW, DZSNSO, ZSNSO);

			}
			
			if (ISNOW < 0) {
				//该调用会更新SH2O的顶层即1层
				COMBINE(Driver.NSnow, Driver.NSoil, ILOC, JLOC,
					ref ISNOW, SH2O, STC, SNICE, SNLIQ,
					DZSNSO, SICE, ref SNOWH, ref SNEQV,
					out PONDING1, out PONDING2);

			}
//			if (ISNOW > 0)
//				throw new Exception("");
			if (ISNOW < 0) {
				DIVIDE(NSNOW, NSOIL, ref ISNOW, STC, SNICE, SNLIQ, DZSNSO);
			}
//			if (ISNOW > 0)
//				throw new Exception("");
			//该函数调用会更新SH2O的顶层，即1层
			SNOWH2O(NSNOW, NSOIL, DT, QSNFRO, QSNSUB, QRAIN, ILOC, JLOC,
				ref ISNOW, DZSNSO, ref SNOWH, ref SNEQV, SNICE,
				SNLIQ, SH2O, SICE, STC,
				out QSNBOT, out PONDING1, out PONDING2);
//			if (double.IsNaN(SNOWH))
//				throw new Exception("");
//			;
			//set empty snow layers to zero
			for (int iz = -NSNOW + 1; iz <= ISNOW; iz++) {  //无雪层以外的层统统设为0
				//do iz = -nsnow+1, isnow
				SNICE[iz] = 0;
				SNLIQ[iz] = 0;
				STC[iz] = 0;
				DZSNSO[iz] = 0;
				ZSNSO[iz] = 0;
			}

			//to obtain equilibrium state of snow in glacier region

			if (SNEQV > 2000) {      // 2000 mm -> maximum water depth
				BDSNOW = SNICE[0] / DZSNSO[0];
				SNOFLOW = (SNEQV - 2000);
				SNICE[0] -= SNOFLOW;
				DZSNSO[0] -= SNOFLOW / BDSNOW;
				SNOFLOW /= DT;
			}

			;   // sum up snow mass for layered snow


			if (ISNOW < 0) {     // MB: only do for multi-layer
				SNEQV = 0;
				for (IZ = ISNOW + 1; IZ <= 0; IZ++)
				//DO IZ = ISNOW+1,0;
					SNEQV += SNICE[IZ] + SNLIQ[IZ];
			}

			;   // Reset ZSNSO and layer thinkness DZSNSO

			for (IZ = ISNOW + 1; IZ <= 0; IZ++)
				DZSNSO[IZ] = -DZSNSO[IZ];

			DZSNSO[1] = ZSOIL[1];
			for (IZ = 2; IZ <= Driver.NSoil; IZ++)
				DZSNSO[IZ] = (ZSOIL[IZ] - ZSOIL[IZ - 1]);
			
			ZSNSO[ISNOW + 1] = DZSNSO[ISNOW + 1];
			
			for (IZ = ISNOW + 2; IZ <= NSOIL; IZ++) {
				ZSNSO[IZ] = ZSNSO[IZ - 1] + DZSNSO[IZ];
			}
			
			for (int iz = ISNOW + 1; iz <= NSOIL; iz++) {
				DZSNSO[iz] *= -1;
			}
		}
		/// <summary>
		/// 该函数已核对过多遍了，没问题了
		/// </summary>
		/// <param name="NSNOW"></param>
		/// <param name="NSOIL"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="ISNOW"></param>
		/// <param name="SH2O"></param>
		/// <param name="STC"></param>
		/// <param name="SNICE"></param>
		/// <param name="SNLIQ"></param>
		/// <param name="DZSNSO"></param>
		/// <param name="SICE"></param>
		/// <param name="SNOWH"></param>
		/// <param name="SNEQV"></param>
		/// <param name="PONDING1"></param>
		/// <param name="PONDING2"></param>
		public static void COMBINE(int NSNOW, int  NSOIL, int ILOC, int JLOC,
			ref int ISNOW, FortDoubleArray SH2O, FortDoubleArray STC, FortDoubleArray SNICE, FortDoubleArray SNLIQ,
			FortDoubleArray DZSNSO, FortDoubleArray SICE, ref double SNOWH, ref double SNEQV,
			out double PONDING1, out double PONDING2)
		{

			// input
			/*
    int, INTENT(IN)      ILOC
    int, INTENT(IN)      JLOC
    int, INTENT(IN)      NSNOW                        //maximum no. of snow layers
    int, INTENT(IN)      NSOIL                        //no. of soil layers

/ input and output

    int,                         INTENT(INOUT)  ISNOW //actual no. of snow layers
    double, DIMENSION(       1:NSOIL), INTENT(INOUT)  SH2O  //soil liquid moisture (m3/m3)
    double, DIMENSION(       1:NSOIL), INTENT(INOUT)  SICE  //soil ice moisture (m3/m3)
    double, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT)  STC   //snow layer temperature [k]
    double, DIMENSION(-NSNOW+1:    0), INTENT(INOUT)  SNICE //snow layer ice [mm]
    double, DIMENSION(-NSNOW+1:    0), INTENT(INOUT)  SNLIQ //snow layer liquid water [mm]
    double, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT)  DZSNSO//snow layer depth [m]
    double,                            INTENT(INOUT)  sneqv //snow water equivalent [m]
    double,                            INTENT(INOUT)  snowh //snow depth [m]
    double,                            INTENT(OUT)  PONDING1
    double,                            INTENT(OUT)  PONDING2
			 */
			// local variables:

			int I, J, K, L;               // node indices
			
			int ISNOW_OLD = 0;           // number of top snow layer
			int MSSI;                 // node index
			int NEIBOR;               // adjacent node selected for combination
			double ZWICE;               // total ice mass in snow
			double ZWLIQ;               // total liquid water in snow
			// minimum of top snow layer
			//    DATA DZMIN /0.045, 0.05, 0.2/
			FortDoubleArray DZMIN = new FortDoubleArray(1, 3);
			DZMIN.data = new double[]{ 0.025, 0.025, 0.1 };  // MB: change limit
			ISNOW_OLD = ISNOW;
			PONDING1 = 0;
			PONDING2 = 0;

			//DO J = ISNOW_OLD+1,0
			for (J = ISNOW_OLD + 1; J <= 0; J++) {
				if (SNICE[J] <= 0.1) {
					if (J != 0) {
						SNLIQ[J + 1] += SNLIQ[J];
						SNICE[J + 1] += SNICE[J];
					} else {
						if (ISNOW_OLD < -1) {     // MB/KM: change to ISNOW
							SNLIQ[J - 1] += SNLIQ[J];
							SNICE[J - 1] += SNICE[J];
						} else {
							if (SNICE[J] >= 0) {
								PONDING1 = SNLIQ[J];    // ISNOW WILL GET SET TO ZERO BELOW; PONDING1 WILL GET
								SNEQV = SNICE[J];       // ADDED TO PONDING FROM PHASECHANGE PONDING SHOULD BE
								SNOWH = DZSNSO[J];      // ZERO HERE BECAUSE IT WAS CALCULATED FOR THIN SNOW
							} else {   // SNICE OVER-SUBLIMATED EARLIER
								PONDING1 = SNLIQ[J] + SNICE[J];
								if (PONDING1 < 0) {   // if SNICE AND SNLIQ SUBLIMATES REMOVE FROM SOIL
									SICE[1] = Math.Max(0.0, SICE[1] + PONDING1 / (DZSNSO[1] * 1000));
									PONDING1 = 0.0;
								}
								SNEQV = 0.0;
								SNOWH = 0.0;
							}
							SNLIQ[J] = 0.0;
							SNICE[J] = 0.0;
							DZSNSO[J] = 0.0;
						}
						//                SH2O[1] = SH2O[1]+SNLIQ[J]/(DZSNSO[1]*1000)
						//                SICE[1] = SICE[1]+SNICE[J]/(DZSNSO[1]*1000)
					}
					
					// shift all elements above this down by one.
					if (J > ISNOW + 1 && ISNOW < -1) {
						for (I = J; I >= ISNOW + 2; I--) {
							//DO I = J, ISNOW+2, -1
							STC[I] = STC[I - 1];
							SNLIQ[I] = SNLIQ[I - 1];
							SNICE[I] = SNICE[I - 1];
							DZSNSO[I] = DZSNSO[I - 1];
						}
					}
					ISNOW += 1;
				} //}
			} //end for

			// to conserve water in case of too large surface sublimation
			//土壤部分从0开始
			if (SICE[1] < 0) {
				SH2O[1] += SICE[1];
				SICE[1] = 0;
				if (SH2O[1] < 0)
					throw new Exception();
			}

			if (ISNOW == 0) {
				return;  // MB: get out if no longer multi-layer
			}

			SNEQV = 0;
			SNOWH = 0;
			ZWICE = 0;
			ZWLIQ = 0;

			for (J = ISNOW + 1; J <= 0; J++) {
				//DO J = ISNOW+1,0
				SNEQV += SNICE[J] + SNLIQ[J];
				SNOWH += DZSNSO[J];
				ZWICE += SNICE[J];
				ZWLIQ += SNLIQ[J];
			}
//			if (double.IsNaN(SNOWH)) {
//				throw new Exception();
//			}
			// check the snow depth - all snow gone
			// the liquid water assumes ponding on soil surface.

			if (SNOWH < 0.025 && ISNOW < 0) {  // MB: change limit
				//       if (SNOWH < 0.05 && ISNOW < 0 )
				ISNOW = 0;
				SNEQV = ZWICE;
				PONDING2 = ZWLIQ;           // LIMIT OF ISNOW < 0 MEANS INPUT PONDING
				if (SNEQV <= 0)
					SNOWH = 0; // SHOULD BE ZERO; SEE ABOVE
			}
			

			//       if (SNOWH < 0.05 )
			//          ISNOW  = 0
			//          SNEQV = ZWICE
			//          SH2O[1] = SH2O[1] + ZWLIQ / (DZSNSO[1] * 1000)
			//          if(SNEQV <= 0) SNOWH = 0.
			//       }

			// check the snow depth - snow layers combined
			if (ISNOW < -1) {

				ISNOW_OLD = ISNOW;
				MSSI = 1;

				//DO I = ISNOW_OLD+1,0
				for (I = ISNOW_OLD + 1; I <= 0; I++) {
					if (DZSNSO[I] < DZMIN[MSSI]) {

						if (I == ISNOW + 1) {
							NEIBOR = I + 1;
						} else if (I == 0) {
							NEIBOR = I - 1;
						} else {
							NEIBOR = I + 1;
							if ((DZSNSO[I - 1] + DZSNSO[I]) < (DZSNSO[I + 1] + DZSNSO[I]))
								NEIBOR = I - 1;
						}

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
						COMBO(ref dzsnso, ref snliq, ref snice, ref stc, DZSNSO[L], SNLIQ[L], SNICE[L], STC[L]);
//						if (double.IsNaN(dzsnso))
//							throw new Exception("");
						DZSNSO[J] = dzsnso;
						SNLIQ[J] = snliq;
						SNICE[J] = snice;
						STC[J] = stc;
						// Now shift all elements above this down one.
						if (J - 1 > ISNOW + 1) {
							//DO K = J-1, ISNOW+2, -1
							for (K = J - 1; K >= ISNOW + 2; K--) {
								STC[K] = STC[K - 1];
								SNICE[K] = SNICE[K - 1];
								SNLIQ[K] = SNLIQ[K - 1];
								DZSNSO[K] = DZSNSO[K - 1];
							}
						}
						// Decrease the number of snow layers
						ISNOW += 1;
						if (ISNOW >= -1)
							break;
					} else {

						// The layer thickness is greater than the prescribed minimum value
						MSSI += 1;

					}
				}
			}
			
		}
		// DIVIDE
		/// <summary>
		/// 核查过了多遍，没问题。2022年11月16日
		/// </summary>
		/// <param name="NSNOW"></param>
		/// <param name="NSOIL"></param>
		/// <param name="DT"></param>
		/// <param name="QSNFRO"></param>
		/// <param name="QSNSUB"></param>
		/// <param name="QRAIN"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="ISNOW"></param>
		/// <param name="DZSNSO"></param>
		/// <param name="SNOWH"></param>
		/// <param name="SNEQV"></param>
		/// <param name="SNICE"></param>
		/// <param name="SNLIQ"></param>
		/// <param name="SH2O"></param>
		/// <param name="SICE"></param>
		/// <param name="STC"></param>
		/// <param name="QSNBOT"></param>
		/// <param name="PONDING1"></param>
		/// <param name="PONDING2"></param>
		static void SNOWH2O(int NSNOW, int NSOIL, double DT, double QSNFRO, double QSNSUB,  //in
			double     QRAIN, int ILOC, int JLOC,                  //in
			ref int    ISNOW, FortDoubleArray DZSNSO, ref double SNOWH, ref double SNEQV, FortDoubleArray SNICE,  //inout
			FortDoubleArray   SNLIQ, FortDoubleArray SH2O, FortDoubleArray SICE, FortDoubleArray STC,          //inout
			out double   QSNBOT, out double PONDING1, out double PONDING2)          //out
		{
			// ----------------------------------------------------------------------
			// Renew the mass of ice lens (SNICE) and liquid (SNLIQ) of the
			// surface snow layer resulting from sublimation (frost) / evaporation (dew)

			// input

			//   INTEGER,                         INTENT(IN)    :: ILOC   //grid index
			//   INTEGER,                         INTENT(IN)    :: JLOC   //grid index
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
//	  REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O   //soil liquid moisture (m3/m3)
			//   REAL, DIMENSION(       1:NSOIL), INTENT(INOUT) :: SICE   //soil ice moisture (m3/m3)
			//   REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC    //snow layer temperature [k]
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
			//   REAL :: PONDING1, PONDING2
			// ----------------------------------------------------------------------

			//for the case when SNEQV becomes '0' after 'COMBINE'
			PONDING1 = 0;
			PONDING2 = 0;
			if (SNEQV == 0) {

				SICE[1] += (QSNFRO - QSNSUB) * DT / (DZSNSO[1] * 1000);  // Barlage: SH2O->SICE v3.6
				if (SICE[1] < 0) {
					SH2O[1] += SICE[1];
					SICE[1] = 0;					
				}
			}

			// for shallow snow without a layer
			// snow surface sublimation may be larger than existing snow mass. To conserve water,
			// excessive sublimation is used to reduce soil water. Smaller time steps would tend
			// to aviod this problem.

			if (ISNOW == 0 && SNEQV > 0) {
				double TEMP = SNEQV;
				SNEQV += -QSNSUB * DT + QSNFRO * DT;
				double PROPOR = SNEQV / TEMP;
				SNOWH = Math.Max(0, PROPOR * SNOWH);

				if (SNEQV < 0) {
					SICE[1] += SNEQV / (DZSNSO[1] * 1000);
					SNEQV = 0;
					SNOWH = 0;
				}
				if (SICE[1] < 0) {
					SH2O[1] += SICE[1];
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
					COMBINE(NSNOW, NSOIL, ILOC, JLOC, //in
						ref ISNOW, SH2O, STC, SNICE, SNLIQ,  //inout
						DZSNSO, SICE, ref SNOWH, ref SNEQV,          //inout
						out PONDING1, out PONDING2);                      //out
				}
				//KWM:  Subroutine COMBINE can change ISNOW to make it 0 again?
				if (ISNOW < 0) { //KWM added this IF statement to prevent out-of-bounds array references
					SNLIQ[ISNOW + 1] += QRAIN * DT;
					SNLIQ[ISNOW + 1] = Math.Max(0, SNLIQ[ISNOW + 1]);
				}

			} //KWM  -- Can the } be moved toward the end of the subroutine (Just set QSNBOT=0)?

			// Porosity and partial volume

			//KWM Looks to me like loop index / IF test can be simplified.

			//DO J = -NSNOW+1, 0
			FortDoubleArray VOL_ICE = new FortDoubleArray(1 - NSNOW, 0);
			FortDoubleArray VOL_LIQ = new FortDoubleArray(1 - NSNOW, 0);
			FortDoubleArray EPORE = new FortDoubleArray(1 - NSNOW, 0);
			for (int J = -NSNOW + 1; J <= 0; J++) {
				if (J >= ISNOW + 1) {
					VOL_ICE[J] = Math.Min(1, SNICE[J] / (DZSNSO[J] * NoahMP.DENICE));
					EPORE[J] = 1 - VOL_ICE[J];
					VOL_LIQ[J] = Math.Min(EPORE[J], SNLIQ[J] / (DZSNSO[J] * NoahMP.DENH2O));
				}
			}

			double QIN = 0;
			double QOUT = 0;

			//KWM Looks to me like loop index / IF test can be simplified.

			//DO J = -NSNOW+1, 0
			for (int J = -NSNOW + 1; J <= 0; J++) {
				if (J >= ISNOW + 1) {
					SNLIQ[J] += QIN;
					if (J <= -1) {
						if (EPORE[J] < 0.05 || EPORE[J + 1] < 0.05) {
							QOUT = 0;
						} else {
							QOUT = Math.Max(0, (VOL_LIQ[J] - NoahMP.SSI * EPORE[J]) * DZSNSO[J]);
							QOUT = Math.Min(QOUT, (1 - VOL_ICE[J + 1] - VOL_LIQ[J + 1]) * DZSNSO[J + 1]);
						}
					} else {
						QOUT = Math.Max(0, (VOL_LIQ[J] - NoahMP.SSI * EPORE[J]) * DZSNSO[J]);
					}
					QOUT *= 1000;
					SNLIQ[J] -= QOUT;
					QIN = QOUT;
				}
			}

			// Liquid water from snow bottom to soil
			
			QSNBOT = QOUT / DT;           // mm/s

		}
		/// <summary>
		/// melting/freezing of snow water and soil water
		/// 2022年11月16日又检查过
		/// </summary>
		/// <param name="NSNOW"></param>
		/// <param name="NSOIL"></param>
		/// <param name="ISNOW"></param>
		/// <param name="DT"></param>
		/// <param name="FACT"></param>
		/// <param name="DZSNSO"></param>
		/// <param name="HCPCT"></param>
		/// <param name="IST"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="STC"></param>
		/// <param name="SNICE"></param>
		/// <param name="SNLIQ"></param>
		/// <param name="SNEQV"></param>
		/// <param name="SNOWH"></param>
		/// <param name="SMC"></param>
		/// <param name="SH2O"></param>
		/// <param name="QMELT"></param>
		/// <param name="IMELT"></param>
		/// <param name="PONDING"></param>
		static void	PHASECHANGE(GridCell cell, int NSNOW, int NSOIL, int ISNOW, double DT, FortDoubleArray FACT, // !in
			FortDoubleArray DZSNSO, FortDoubleArray HCPCT, int IST, int ILOC, int JLOC, // !in
			FortDoubleArray STC, FortDoubleArray SNICE, FortDoubleArray SNLIQ, ref double SNEQV, ref double SNOWH, // !inout
			FortDoubleArray SMC, FortDoubleArray SH2O,                            // !inout
			out double QMELT, FortIntArray IMELT, out double PONDING)                    // !out
		{

			// inputs

			//  INTEGER, INTENT(IN)                             :: ILOC   //grid index
			//  INTEGER, INTENT(IN)                             :: JLOC   //grid index
			//  INTEGER, INTENT(IN)                             :: NSNOW  //maximum no. of snow layers [=3]
			//  INTEGER, INTENT(IN)                             :: NSOIL  //No. of soil layers [=4]
			//  INTEGER, INTENT(IN)                             :: ISNOW  //actual no. of snow layers [<=3]
			//  INTEGER, INTENT(IN)                             :: IST    //surface type: 1->soil; 2->lake
			//  REAL, INTENT(IN)                                :: DT     //land model time step (sec)
			//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)     :: FACT   //temporary
			//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)     :: DZSNSO //snow/soil layer thickness [m]
			//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)     :: HCPCT  //heat capacity (J/m3/k)
//
			//// outputs
			//  INTEGER, DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: IMELT  //phase change index
			//  REAL,                               INTENT(OUT) :: QMELT  //snowmelt rate [mm/s]
			//  REAL,                               INTENT(OUT) :: PONDING//snowmelt when snow has no layer [mm]
//
			//// inputs and outputs
//
			//  REAL, INTENT(INOUT) :: SNEQV
//			 REAL, INTENT(INOUT) :: SNOWH
			//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT)  :: STC    //snow/soil layer temperature [k]
			//  REAL, DIMENSION(       1:NSOIL), INTENT(INOUT)  :: SH2O   //soil liquid water [m3/m3]
			//  REAL, DIMENSION(       1:NSOIL), INTENT(INOUT)  :: SMC    //total soil water [m3/m3]
			//  REAL, DIMENSION(-NSNOW+1:0)    , INTENT(INOUT)  :: SNICE  //snow layer ice [mm]
			//  REAL, DIMENSION(-NSNOW+1:0)    , INTENT(INOUT)  :: SNLIQ  //snow layer liquid water [mm]
//
			//// local
//
			//  INTEGER                         :: J         //do loop index
			//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: HM        //energy residual [w/m2]
			//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: XM        //melting or freezing water [kg/m2]
			//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: WMASS0
			//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: WICE0
			//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: WLIQ0
			//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: MICE      //soil/snow ice mass [mm]
			//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: MLIQ      //soil/snow liquid water mass [mm]
			//  REAL, DIMENSION(-NSNOW+1:NSOIL) :: SUPERCOOL //supercooled water in soil (kg/m2)
			//  REAL                            :: HEATR     //energy residual or loss after melting/freezing
			//  REAL                            :: TEMP1     //temporary variables [kg/m2]
			//  REAL                            :: PROPOR
			//  REAL                            :: SMP       //frozen water potential (mm)
			//  REAL                            :: XMF       //total latent heat of phase change

			// ----------------------------------------------------------------------
			// Initialization
//		if (double.IsNaN(STC[0]))
//					throw new Exception("");
			QMELT = 0;
			PONDING = 0;
			double XMF = 0;
			FortDoubleArray SUPERCOOL = new FortDoubleArray(1 - NSNOW, NSOIL);
			FortDoubleArray MICE = new FortDoubleArray(1 - NSNOW, NSOIL);
			FortDoubleArray MLIQ = new FortDoubleArray(1 - NSNOW, NSOIL);
			
			for (int J = ISNOW + 1; J <= 0; J++) {
				MICE[J] = SNICE[J];
				MLIQ[J] = SNLIQ[J];
			}

			//DO J = 1, NSOIL               // soil
			for (int J = 1; J <= NSOIL; J++) {
				MLIQ[J] = SH2O[J] * DZSNSO[J] * 1000;
				MICE[J] = (SMC[J] - SH2O[J]) * DZSNSO[J] * 1000;
			}

			
			FortDoubleArray HM = new FortDoubleArray(1 - NSNOW, NSOIL);
			FortDoubleArray XM = new FortDoubleArray(1 - NSNOW, NSOIL);
			FortDoubleArray WICE0 = new FortDoubleArray(1 - NSNOW, NSOIL);
			FortDoubleArray WLIQ0 = new FortDoubleArray(1 - NSNOW, NSOIL);
			FortDoubleArray WMASS0 = new FortDoubleArray(1 - NSNOW, NSOIL);
			//IMELT=new int[NSOIL+NSNOW];
			//DO J = ISNOW+1,NSOIL       // all layers
			for (int J = ISNOW + 1; J <= NSOIL; J++) {
				IMELT[J] = 0;
				HM[J] = 0;
				XM[J] = 0;
				WICE0[J] = MICE[J];
				WLIQ0[J] = MLIQ[J];
				WMASS0[J] = MICE[J] + MLIQ[J];
			}

			
			if (IST == 1) {
				for (int J = 1; J <= NSOIL; J++) {
					if (NoahMP.OPT_FRZ == 1) {
						if (STC[J] < NoahMP.TFRZ) {
							double SMP = NoahMP.HFUS * (NoahMP.TFRZ - STC[J]) / (NoahMP.GRAV * STC[J]);             //(m)
							SUPERCOOL[J] = cell.SMCMAX * Math.Pow(SMP / cell.PSISAT, -1 / cell.BEXP);
							SUPERCOOL[J] *= DZSNSO[J] * 1000;        //(mm)
						}						
					}
					if (NoahMP.OPT_FRZ == 2) {  //验证完毕，问题不大, 误差来源于STC
						
						SUPERCOOL[J] = FRH2O(cell, STC[J], SMC[J], SH2O[J]);
						SUPERCOOL[J] *= DZSNSO[J] * 1000;        //(mm)
					}
				}
			}
			
			//DO J = ISNOW+1,NSOIL
			for (int J = ISNOW + 1; J <= NSOIL; J++) {
				if (MICE[J] > 0 && STC[J] >= NoahMP.TFRZ) { //melting
					IMELT[J] = 1;
				}
				if (MLIQ[J] > SUPERCOOL[J] && STC[J] < NoahMP.TFRZ) {
					IMELT[J] = 2;
				}

				// if snow exists, but its thickness is not enough to create a layer
				if (ISNOW == 0 && SNEQV > 0 && J == 1) {
					if (STC[J] >= NoahMP.TFRZ) {
						IMELT[J] = 1;
					}
				}
			}
			// Calculate the energy surplus and loss for melting and freezing

			for (int J = ISNOW + 1; J <= NSOIL; J++) {
				if (IMELT[J] > 0) {
					HM[J] = (STC[J] - NoahMP.TFRZ) / FACT[J];
					STC[J] = NoahMP.TFRZ;
				}

				if (IMELT[J] == 1 && HM[J] < 0) {
					HM[J] = 0;
					IMELT[J] = 0;
				}
				if (IMELT[J] == 2 && HM[J] > 0) {
					HM[J] = 0;
					IMELT[J] = 0;
				}
				XM[J] = HM[J] * DT / NoahMP.HFUS;
			}

			// The rate of melting and freezing for snow without a layer, needs more work.

			if (ISNOW == 0 && SNEQV > 0 && XM[1] > 0) {
				double TEMP1 = SNEQV;
				SNEQV = Math.Max(0, TEMP1 - XM[1]);
				double	PROPOR = SNEQV / TEMP1;
				SNOWH = Math.Max(0, PROPOR * SNOWH);
				double HEATR = HM[1] - NoahMP.HFUS * (TEMP1 - SNEQV) / DT;
				if (HEATR > 0) {
					XM[1] = HEATR * DT / NoahMP.HFUS;
					HM[1] = HEATR;
				} else {
					XM[1] = 0;
					HM[1] = 0;
				}
				QMELT = Math.Max(0, (TEMP1 - SNEQV)) / DT;
				XMF = NoahMP.HFUS * QMELT;
				PONDING = TEMP1 - SNEQV;
			}

			//The rate of melting and freezing for snow and soil
			for (int J = ISNOW + 1; J <= NSOIL; J++) {
				if (IMELT[J] > 0 && Math.Abs(HM[J]) > 0) {
					double HEATR = 0;
					if (XM[J] > 0) {
						MICE[J] = Math.Max(0, WICE0[J] - XM[J]);
						HEATR = HM[J] - NoahMP.HFUS * (WICE0[J] - MICE[J]) / DT;
					} else if (XM[J] < 0) {
						if (J <= 0) {                             // snow
							MICE[J] = Math.Min(WMASS0[J], WICE0[J] - XM[J]);
						} else {                                        // soil
							if (WMASS0[J] < SUPERCOOL[J]) {
								MICE[J] = 0;
							} else {
								MICE[J] = Math.Min(WMASS0[J] - SUPERCOOL[J], WICE0[J] - XM[J]);
								MICE[J] = Math.Max(MICE[J], 0.0);
							}
						}
						HEATR = HM[J] - NoahMP.HFUS * (WICE0[J] - MICE[J]) / DT;
						
					}
//					if(double.IsNaN(HEATR))
//						throw new Exception();
					MLIQ[J] = Math.Max(0, WMASS0[J] - MICE[J]);
					if (Math.Abs(HEATR) > 0) {
						STC[J] += FACT[J] * HEATR;
						if (J <= 0) {                             // snow
							if (MLIQ[J] * MICE[J] > 0)
								STC[J] = NoahMP.TFRZ;
						}
					}

					XMF += NoahMP.HFUS * (WICE0[J] - MICE[J]) / DT;

					if (J < 1) {
						QMELT += Math.Max(0, (WICE0[J] - MICE[J])) / DT;
					}
				}
			}//end do

			for (int J = ISNOW + 1; J <= 0; J++) {
				SNLIQ[J] = MLIQ[J];
				SNICE[J] = MICE[J];
			}

			//DO J = 1, NSOIL             // soil
			for (int J = 1; J <= NSOIL; J++) {     //2022年11月17日晚发现了这里缺少等号，解决了连续多日疑惑的一个问题
				SH2O[J] = MLIQ[J] / (1000 * DZSNSO[J]);
				SMC[J] = (MLIQ[J] + MICE[J]) / (1000 * DZSNSO[J]);
			}
		}
		public static void CARBON(int NSNOW, int NSOIL, int VEGTYP, int NROOT, double DT, FortDoubleArray ZSOIL,    //in
			FortDoubleArray   DZSNSO, FortDoubleArray STC, FortDoubleArray SMC, double TV, double TG, double PSN,    //in
			double   FOLN, double SMCMAX, double BTRAN, double APAR, double FVEG, double IGS,   //in
			double  TROOT, double IST, double LAT, double ILOC, double JLOC, int ISURBAN,   //in
			ref double   LFMASS, ref double RTMASS, ref double STMASS, ref double WOOD, ref double STBLCP, ref double FASTCP,  //inout
			out  double  GPP, out double NPP, out double NEE, out double AUTORS, out double HETERS,   //out
			out double TOTSC, out double TOTLB, ref double XLAI, ref double XSAI)                   //out
		{
			
			//初始化下面的输出变量
			TOTSC = 0;
			TOTLB = 0;
			
			if ((VEGTYP == Driver.vegparams.ISWATER)
			    || (VEGTYP == Driver.vegparams.ISBARREN)
			    || (VEGTYP == Driver.vegparams.ISSNOW)
			    || (VEGTYP == Driver.vegparams.ISURBAN)) {
				
				GPP = 0;
				NPP = 0;
				NEE = 0;
				AUTORS = 0;
				HETERS = 0;
				LFMASS = 0;
				RTMASS = 0;
				STMASS = 0;
				WOOD = 0;
				STBLCP = 0;
				FASTCP = 0;
				return;
			}
			//LAPM: leaf area per unit mass [m2/g]，单位质量的叶片面积，与具体植被有关
			double LAPM = Driver.vegparams.SLA[VEGTYP - 1] / 1000;   // m2/kg -> m2/g
			

			//water stress coeficient [-]  (1. for wilting )
			double WSTRES = 1 - BTRAN;
			//root zone soil water [-]
			double WROOT = 0;
			for (int J = 1; J <= NROOT; J++)
				WROOT += SMC[J] / SMCMAX * DZSNSO[J] / (-ZSOIL[NROOT]);
			CO2FLUX_52(NSNOW, NSOIL, VEGTYP, IGS, DT,  //in
				DZSNSO, STC, PSN, TROOT, TV,  //in
				WROOT, WSTRES, FOLN, LAPM,          //in
				LAT, ILOC, JLOC, FVEG,          //in
				ref XLAI, ref  XSAI, ref LFMASS, ref RTMASS, ref STMASS,  //inout
				ref	FASTCP, ref STBLCP, ref WOOD,                  //inout
				out	GPP, out NPP, out NEE, out AUTORS, out HETERS,  //out
				out	TOTSC, out TOTLB);                      //out

		}
		/// <summary>
		/// 刘永和加写的代码
		/// 根据文献：Sibo Zeng, Zaihua Liu, Georg Kaufmann(2019),Sensitivity of the global carbonate weathering carbon-sink flux to 
		/// climate and land-use changes, Nature communications. 
		/// https://www.nature.com/articles/s41467-019-13772-4
		/// </summary>
		/// <param name="NPP">净初级生产力</param>
		/// <param name="RUNOFF">径流量</param>
		/// <param name="TG">地温</param>
		/// <param name="SFCPRS">地表面气压</param>
		/// param name="CCO2">全球二氧化碳浓度，随着排放量逐年变化,2022年时为417ppm</param>
		/// <param name="CCSF">碳酸岩风化碳汇通量</param>
		public static void CarbonSink(double NPP, double RUNOFF, double TG, double SFCPRS, double CCO2, out double CCSF)
		{			
			double pCO2 = CCO2 * 1e-6;//pCO2的单位要求是标准大气压(1atm)。
			double T = TG;
			double logT = Math.Log10(T);
			double TC = TG - 273.16;
			double KH = Math.Pow(10, 108.3865 - 6919.53 / T + 0.01985076 * T - 40.45154 * logT + 669365 / (T * T));
			double K1 = Math.Pow(10, -356.3094 - 0.06091964 * T + 21834.37 / T + 126.8339 * logT - 1684915 / (T * T));
			double K2 = Math.Pow(10, -107.8871 + 5151.79 / T - 0.03252849 * T + 38.92561 * logT - 563713.9 / (T * T));
			double KC = Math.Pow(10, -171.9065 - 0.077993 * T + 2839.319 / T + 71.595 * logT);			
			
			double A = 0.4883 + TC * 8.074E-4;
			double B = 0.3241 + TC * 1.6E-4; 
			double A2 = 1.03e6;
			double pCO2_SOIL = pCO2 + A2 * 0.75 * NPP / (TG * TG); //Eq.(3),土壤CO2分压		
			//**********************求取Ca开始**************************
			double Ca = 1e-4;  //任意设置的一个Ca离子初始浓度，从此值开始迭代，网上查到的结果是1.5e-4mol/L。碳酸钙的溶解度为0.0014(g/100g水)，即一升水有0.014g/L,除以分子量100后，得离子浓度为0.00014mol/L
			for (int i = 0; i < 6; i++) {
				double IonicStrength = Ca * 3; //这一个参数现在还不确定取几
				double rootI = Math.Sqrt(IonicStrength); 
				//下面使用的是Debye-Huckel equation
				double gammaCa = Math.Pow(10, (-A * 4 * rootI / (1 + B * 5.0 * rootI) + 0.165 * IonicStrength));   //mol/m3
				double gammaHC03 = Math.Pow(10, (-A * rootI / (1 + B * 5.4 * rootI)));
				double gammaOH = Math.Pow(10, (-A * rootI / (1 + B * 3.5 * rootI)));
				double gammaH = Math.Pow(10, (-A * rootI / (1 + B * 9 * rootI)));
				double gammaC03 = Math.Pow(10, (-A * rootI / (1 + B * 5.4 * rootI)));
				Ca = Math.Pow(K1 * KC * KH * pCO2_SOIL / (4 * K2 * gammaCa * gammaHC03), 1.0 / 3); //MOL/m3
			}
			//**********************求取Ca结束**************************
			double CaEq = Ca * 1000;  //unit: mol/m3    计算钙离子浓度,据郭芳，本人换算的该值为0.002mol/L
			//CaEq2是采用近似公式
			double CaEq2 = 10.75 * (1 - 0.0139 * TC) * Math.Pow(pCO2_SOIL, 1.0 / 3); //unit: 微mol/ML，或mmol/L			
			CCSF = 12 * RUNOFF * CaEq;			
		}
		
		
		
//		public static void RWEQ(double SC,double CL,double OM,double L,double dH,double Crr,double z,double El,double T,double SR,double DT,double Prosnow,double EF)
//		{
//			double COG=Math.Exp(-0.0483*SC);
//			double SCF=1.0/(1+0.0066*CL*CL+0.021*OM*OM);
//			double Kr=0.2*dH*dH/L;
//			double Kp=Math.Exp(1.86*Kr-2.41*Math.Pow(Kr,0.934)-0.124*Crr);		
//			
//			double rho=348.01*(1.013-0.1183*El+0.0048*El*El)/T;  //要求T是绝对温度
//			
//			double ETp=0.0162*SR/58.5*(DT+17.8);
//			double SW=ETp-R
//			double SD=1-Prosnow;
//			double WF=0;
//			//EF为土壤可蚀性因子
//			double Qmax=109.8*WF*EF*SCF*Kp*COG;
//			double S=150.71*Math.Exp(WF*EF*SCF*Kp*COG,-0.3711);
//			double SL=2*z/(S*S)*Qmax*Math.Exp(-(z/S)*(z/S));
//			
//		}
		/// <summary>
		/// 通过单位线计算坡面汇流量
		/// </summary>
		/// <param name="runoff"></param>
		public static void CalHillSlopeFlow(GridCell cell, double runoff)
		{
			switch (NoahMP.OPT_UHG) {
				case 0:
					{
						//相当于什么也不做
						cell.Q_HillSlope = 0;
						break;
					}
				case 1:
					{
						cell.Q_HillSlope = runoff;
						break;
					}
				case 2:
					{
						for (int i = 1; i < 10; i++) {
							cell.HistoricalRunoff[i] = cell.HistoricalRunoff[i - 1];
						}
						cell.HistoricalRunoff[0] = runoff;
						double sum = 0;
						double sumUnit = 0;
						for (int tau = 1; tau < 10; tau++) {
							double unit = cell.NashUnit(tau);
							sumUnit += unit;
							sum += cell.HistoricalRunoff[tau] * unit;
						}
						sum += cell.HistoricalRunoff[0] * (1 - sumUnit); //这句分开写是为了使单位权重之和为1.
						cell.Q_HillSlope = sum;
						break;
					}
				case 3:
					{
						cell.HistoricalRunoff[0] += runoff;
						double Q = cell.HistoricalRunoff[0] * 0.5; //系数为0.5
						cell.Q_HillSlope = Q;
						cell.HistoricalRunoff[0] = cell.HistoricalRunoff[0] - Q;
						break;
					}
			}
		}

		/// 已检查过了，应该无错了,这是之前写的。后来在2021年9月时又检查过了一遍
		/// </summary>
		/// <param name="NSNOW"></param>
		/// <param name="NSOIL"></param>
		/// <param name="DT"></param>
		/// <param name="SICE"></param>
		/// <param name="ZSOIL"></param>
		/// <param name="STC"></param>
		/// <param name="WCND"></param>
		/// <param name="FCRMAX"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="SH2O"></param>
		/// <param name="ZWT"></param>
		/// <param name="WA"></param>
		/// <param name="WT"></param>
		/// <param name="QIN"></param>
		/// <param name="QDIS"></param>
		static void GROUNDWATER(GridCell cell, int NSNOW, int NSOIL, double DT, FortDoubleArray SICE, FortDoubleArray ZSOIL, // !in
			FortDoubleArray STC, FortDoubleArray WCND, double FCRMAX, int ILOC, int JLOC,  //in
			FortDoubleArray  SH2O, ref double ZWT, ref double WA, ref double WT,          //inout
			out double  QIN, out double QDIS)                           //out
		{

			//  INTEGER,                         INTENT(IN) :: ILOC  !grid index
			//  INTEGER,                         INTENT(IN) :: JLOC  !grid index
			//  INTEGER,                         INTENT(IN) :: NSNOW !maximum no. of snow layers
			//  INTEGER,                         INTENT(IN) :: NSOIL !no. of soil layers
			//  REAL,                            INTENT(IN) :: DT    !timestep [sec]
			//  REAL,                            INTENT(IN) :: FCRMAX!maximum FCR (-)
			//  REAL, DIMENSION(       1:NSOIL), INTENT(IN) :: SICE  !soil ice content [m3/m3]
			//  REAL, DIMENSION(       1:NSOIL), INTENT(IN) :: ZSOIL !depth of soil layer-bottom [m]
			//  REAL, DIMENSION(       1:NSOIL), INTENT(IN) :: WCND  !hydraulic conductivity (m/s)
			//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: STC   !snow/soil temperature (k)
//
			//! input and output
			//  REAL, DIMENSION(    1:NSOIL), INTENT(INOUT) :: SH2O  !liquid soil water [m3/m3]
			//  REAL,                         INTENT(INOUT) :: ZWT   !the depth to water table [m]
			//  REAL,                         INTENT(INOUT) :: WA    !water storage in aquifer [mm]
			//  REAL,                         INTENT(INOUT) :: WT    !water storage in aquifer
			//                                                           !+ saturated soil [mm]
			//! output
			//  REAL,                           INTENT(OUT) :: QIN   !groundwater recharge [mm/s]
			//  REAL,                           INTENT(OUT) :: QDIS  !groundwater discharge [mm/s]
//
			//! local
			//  REAL                                        :: FFF   !runoff decay factor (m-1)
			//  REAL                                        :: RSBMX !baseflow coefficient [mm/s]
			//  INTEGER                                     :: IZ    !do-loop index
			//  INTEGER                                     :: IWT   !layer index above water table layer
//			  REAL,  DIMENSION(    1:NSOIL)               :: DZMM  !layer thickness [mm]
			//  REAL,  DIMENSION(    1:NSOIL)               :: ZNODE !node depth [m]
			//  REAL,  DIMENSION(    1:NSOIL)               :: MLIQ  !liquid water mass [kg/m2 or mm]
			//  REAL,  DIMENSION(    1:NSOIL)               :: EPORE !effective porosity [-]
			//  REAL,  DIMENSION(    1:NSOIL)               :: HK    !hydraulic conductivity [mm/s]
			//  REAL,  DIMENSION(    1:NSOIL)               :: SMC   !total soil water  content [m3/m3]
			//  REAL(KIND=8)                                :: S_NODE!degree of saturation of IWT layer
			//  REAL                                        :: DZSUM !cumulative depth above water table [m]
			//  REAL                                        :: SMPFZ !matric potential (frozen effects) [mm]
			//  REAL                                        :: KA    !aquifer hydraulic conductivity [mm/s]
			//  REAL                                        :: WH_ZWT!water head at water table [mm]
			//  REAL                                        :: WH    !water head at layer above ZWT [mm]
			//  REAL                                        :: WS    !water used to fill air pore [mm]
			//  REAL                                        :: WTSUB !sum of HK*DZMM
			//  REAL                                        :: WATMIN!minimum soil vol soil moisture [m3/m3]
			//  REAL                                        :: XS    !excessive water above saturation [mm]
			double ROUS = 0.2;    //!specific yield [-]
			
			//                                                               !0.0-close to free drainage
			// -------------------------------------------------------------
			QDIS = 0.0;
			QIN = 0.0;
			FortDoubleArray DZMM = new FortDoubleArray(1, Driver.NSoil);
			FortDoubleArray ZNODE = new FortDoubleArray(1, Driver.NSoil);
			FortDoubleArray MLIQ = new FortDoubleArray(1, Driver.NSoil);
			FortDoubleArray HK = new FortDoubleArray(1, Driver.NSoil);
			FortDoubleArray EPORE = new FortDoubleArray(1, Driver.NSoil);
			// Derive layer-bottom depth in [mm]
			//KWM:  Derive layer thickness in mm

			DZMM[1] = -ZSOIL[1] * 1000;
			for (int IZ = 2; IZ <= NSOIL; IZ++) {
				DZMM[IZ] = 1000 * (ZSOIL[IZ - 1] - ZSOIL[IZ]);
			}

			// Derive node (middle) depth in [m]
			//KWM:  Positive number, depth below ground surface in m
			ZNODE[1] = -ZSOIL[1] / 2;
			for (int IZ = 2; IZ <= NSOIL; IZ++) {
				ZNODE[IZ] = -ZSOIL[IZ - 1] + 0.5 * (ZSOIL[IZ - 1] - ZSOIL[IZ]);
			}

			// Convert volumetric soil moisture "sh2o" to mass

			FortDoubleArray SMC = new FortDoubleArray(1, NSOIL);
			for (int IZ = 1; IZ <= NSOIL; IZ++) {
				SMC[IZ] = SH2O[IZ] + SICE[IZ];
				MLIQ[IZ] = SH2O[IZ] * DZMM[IZ];
				EPORE[IZ] = Math.Max(0.01, cell.SMCMAX - SICE[IZ]);
				HK[IZ] = 1E3 * WCND[IZ];
			}

			// The layer index of the first unsaturated layer,
			// i.e., the layer right above the water table

			int IWT = NSOIL;
			for (int IZ = 2; IZ <= NSOIL; IZ++) {
				if (ZWT <= -ZSOIL[IZ]) {
					IWT = IZ - 1;
					break;
				}
			}

			// Groundwater discharge [mm/s]

			double FFF = 6.0;
			double RSBMX = 5.0;
			QDIS = (1.0 - FCRMAX) * RSBMX * Math.Exp(-NoahMP.TIMEAN) * Math.Exp(-FFF * (ZWT - 2.0));  //向外排出量,由于临时水位ZWT决定
//			Console.WriteLine("FCRMAX="+FCRMAX.ToString("E")+"TIMEAN="+NoahMP.TIMEAN+"	ZWT="+ZWT.ToString("E"));

			// Matric potential at the layer above the water table

			double S_NODE = Math.Min(1.0, SMC[IWT] / cell.SMCMAX);  //唯一由SMC影响QIN、ZWT的语句
			S_NODE = Math.Max(S_NODE, 0.01);
			double SMPFZ = -cell.PSISAT * 1000 * Math.Pow(S_NODE, -cell.BEXP);  // m --> mm
			double CMIC = 0.20;   //!microprore content (0.0-1.0)
			SMPFZ = Math.Max(-120000.0, CMIC * SMPFZ);

			// Recharge rate Qin to groundwater

			double KA = HK[IWT];
			double WH_ZWT = -ZWT * 1E3;                          //(mm)
			double WH = SMPFZ - ZNODE[IWT] * 1E3;             //(mm)
			QIN = -KA * (WH_ZWT - WH) / ((ZWT - ZNODE[IWT]) * 1E3);  //补给量
			QIN = Math.Max(-10.0 / DT, Math.Min(10 / DT, QIN));

			// Water storage in the aquifer + saturated soil

			WT += (QIN - QDIS) * DT;     //!(mm)  更新地下水储水量
			if (IWT == NSOIL) {
				WA += (QIN - QDIS) * DT;    //!(mm)
				WT = WA;
				ZWT = (-ZSOIL[NSOIL] + 25) - WA / 1000 / ROUS;      //!(m)
				MLIQ[NSOIL] -= QIN * DT;       // [mm]  更新SH2O，最顶层减去补给量
				MLIQ[NSOIL] += Math.Max(0, (WA - 5000));
				WA = Math.Min(WA, 5000);
			} else {
				if (IWT == NSOIL - 1) {
					ZWT = -ZSOIL[NSOIL] - (WT - ROUS * 1000 * 25) / (EPORE[NSOIL]) / 1000;
				} else {
					double	WS = 0;   // water used to fill soil air pores
					//  DO IZ = IWT+2,NSOIL
					for (int IZ = IWT + 2; IZ <= NSOIL; IZ++) {
						WS += EPORE[IZ] * DZMM[IZ];
					}
					ZWT = -ZSOIL[IWT + 1] - (WT - ROUS * 1000 * 25 - WS) / EPORE[IWT + 1] / 1000;
				}

				double WTSUB = 0;
				for (int IZ = 1; IZ <= NSOIL; IZ++) {
					WTSUB += HK[IZ] * DZMM[IZ];
				}

				// Removing subsurface runoff
				for (int IZ = 1; IZ <= NSOIL; IZ++) {
					MLIQ[IZ] -= QDIS * DT * HK[IZ] * DZMM[IZ] / WTSUB;
				}
			}
			
			
			ZWT = Math.Max(1.5, ZWT);
//
			// Limit MLIQ to be greater than or equal to watmin.
			// Get water needed to bring MLIQ equal WATMIN from lower layer.
//
			double WATMIN = 0.01;
			//DO IZ = 1, NSOIL-1
			double XS = 0;
			for (int IZ = 1; IZ <= NSOIL - 1; IZ++) {
				if (MLIQ[IZ] < 0) {
					XS = WATMIN - MLIQ[IZ];
				} else {
					XS = 0;
				}
				MLIQ[IZ] += XS;
				MLIQ[IZ + 1] -= XS;
			}

			int iz = NSOIL;
			if (MLIQ[iz] < WATMIN)
				XS = WATMIN - MLIQ[iz];
			else
				XS = 0;
			
			MLIQ[iz] += XS;
			WA -= XS;
			WT -= XS;
			for (int IZ = 1; IZ <= NSOIL; IZ++) {
				SH2O[IZ] = MLIQ[IZ] / DZMM[IZ];
			}
			
		}
		/// <summary>
		/// NoahMP3中暂不提供有关这个函数的选项
		/// 核查过一轮,后又检查过一轮
		/// </summary>
		/// <param name="NSNOW"></param>
		/// <param name="NSOIL"></param>
		/// <param name="ZSOIL"></param>
		/// <param name="DT"></param>
		/// <param name="DZSNSO"></param>
		/// <param name="SMCEQ"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="SMC"></param>
		/// <param name="WTD"></param>
		/// <param name="SMCWTD"></param>
		/// <param name="RECH"></param>
		/// <param name="QDRAIN"></param>
		static void SHALLOWWATERTABLE(GridCell cell, int NSNOW, int NSOIL, FortDoubleArray ZSOIL, double DT,  //in
			FortDoubleArray  DZSNSO, FortDoubleArray SMCEQ, int ILOC, int JLOC,  //in
			FortDoubleArray  SMC, ref double WTD, ref double SMCWTD, out double RECH, ref double QDRAIN)  //inout
		{
			// ----------------------------------------------------------------------
			//Diagnoses water table depth and computes recharge when the water table is within the resolved soil layers,
			//according to the Miguez-MachoFan scheme
			// input
			//  INTEGER,                         INTENT(IN) :: NSNOW //maximum no. of snow layers
			//  INTEGER,                         INTENT(IN) :: NSOIL //no. of soil layers
			//  INTEGER,                         INTENT(IN) :: ILOC,JLOC
			//  REAL,                            INTENT(IN) :: DT
			//  REAL, DIMENSION(       1:NSOIL), INTENT(IN) :: ZSOIL //depth of soil layer-bottom [m]
			//  REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DZSNSO // snow/soil layer thickness [m]
			//  REAL,  DIMENSION(      1:NSOIL), INTENT(IN) :: SMCEQ  //equilibrium soil water  content [m3/m3]
//
			//// input and output
			//  REAL,  DIMENSION(      1:NSOIL), INTENT(INOUT) :: SMC   //total soil water  content [m3/m3]
			//  REAL,                         INTENT(INOUT) :: WTD   //the depth to water table [m]
			//  REAL,                         INTENT(INOUT) :: SMCWTD   //soil moisture between bottom of the soil and the water table [m3/m3]
			//  REAL,                         INTENT(OUT) :: RECH // groundwater recharge (net vertical flux across the water table), positive up
			//  REAL,                         INTENT(INOUT) :: QDRAIN
//
			//// local
			//  INTEGER                                     :: IZ    //do-loop index
			//  INTEGER                                     :: IWTD   //layer index above water table layer
			//  INTEGER                                     :: KWTD   //layer index where the water table layer is
			//  REAL                                        :: WTDOLD
			//  REAL                                        :: DZUP
			//  REAL                                        :: SMCEQDEEP
			double SMCMAX = cell.SMCMAX;
			RECH = 0;
			FortDoubleArray ZSOIL0 = new FortDoubleArray(0, NSOIL);// new double[NSOIL];
			for (int i = 1; i <= NSOIL; i++) {
				ZSOIL0[i] = ZSOIL[i];
			}
			
			ZSOIL0[0] = 0;

			//find the layer where the water table is
			//     DO IZ=NSOIL,1,-1
			int IZ = 0;
			for (IZ = NSOIL; IZ >= 1; IZ--) {
				if (WTD + 1E-6 < ZSOIL0[IZ])
					break;
			}
			int IWTD = IZ;

			int KWTD = IWTD + 1;  //layer where the water table is
			if (KWTD <= NSOIL) {   //wtd in the resolved layers
				double	WTDOLD = WTD;
				if (SMC[KWTD] > SMCEQ[KWTD]) {

					if (SMC[KWTD] == SMCMAX) { //wtd went to the layer above
						WTD = ZSOIL0[IWTD];
						RECH = -(WTDOLD - WTD) * (SMCMAX - SMCEQ[KWTD]);
						IWTD -= 1;
						KWTD -= 1;
						if (KWTD >= 1) {
							if (SMC[KWTD] > SMCEQ[KWTD]) {
								WTDOLD = WTD;
								WTD = Math.Min((SMC[KWTD] * DZSNSO[KWTD]
								- SMCEQ[KWTD] * ZSOIL0[IWTD] + SMCMAX * ZSOIL0[KWTD]) /
								(SMCMAX - SMCEQ[KWTD]), ZSOIL0[IWTD]);
								RECH -= (WTDOLD - WTD) * (SMCMAX - SMCEQ[KWTD]);
							}
						}
					} else { //wtd stays in the layer
						WTD = Math.Min((SMC[KWTD] * DZSNSO[KWTD] - SMCEQ[KWTD] * ZSOIL0[IWTD] + SMCMAX * ZSOIL0[KWTD]) / (SMCMAX - SMCEQ[KWTD]), ZSOIL0[IWTD]);
						RECH = -(WTDOLD - WTD) * (SMCMAX - SMCEQ[KWTD]);
					}
				} else {    //wtd has gone down to the layer below
					WTD = ZSOIL0[KWTD];
					RECH = -(WTDOLD - WTD) * (SMCMAX - SMCEQ[KWTD]);
					KWTD += 1;
					IWTD += 1;
					//wtd crossed to the layer below. Now adjust it there
					if (KWTD <= NSOIL) {
						WTDOLD = WTD;
						if (SMC[KWTD] > SMCEQ[KWTD]) {
							WTD = Math.Min((SMC[KWTD] * DZSNSO[KWTD]
							- SMCEQ[KWTD] * ZSOIL0[IWTD] + SMCMAX * ZSOIL0[KWTD]) /
							(SMCMAX - SMCEQ[KWTD]), ZSOIL0[IWTD]);
						} else {
							WTD = ZSOIL0[KWTD];
						}
						RECH -= (WTDOLD - WTD) * (SMCMAX - SMCEQ[KWTD]);
					} else {
						WTDOLD = WTD;
						//restore smoi to equilibrium value with water from the ficticious layer below
						//                   SMCWTD=SMCWTD-(SMCEQ[NSOIL]-SMC[NSOIL])
						//                   QDRAIN = QDRAIN - 1000 * (SMCEQ[NSOIL]-SMC[NSOIL]) * DZSNSO[NSOIL] / DT
						//                   SMC[NSOIL]=SMCEQ[NSOIL]
						//adjust wtd in the ficticious layer below
						double SMCEQDEEP = SMCMAX * Math.Pow(-cell.PSISAT / (-cell.PSISAT - DZSNSO[NSOIL]), 1 / cell.BEXP);
						WTD = Math.Min((SMCWTD * DZSNSO[NSOIL]
						- SMCEQDEEP * ZSOIL0[NSOIL] + SMCMAX * (ZSOIL0[NSOIL] - DZSNSO[NSOIL])) /
						(SMCMAX - SMCEQDEEP), ZSOIL0[NSOIL]);
						RECH -= (WTDOLD - WTD) * (SMCMAX - SMCEQDEEP);
					}
				}
			} else if (WTD >= ZSOIL0[NSOIL] - DZSNSO[NSOIL]) {
				//if wtd was already below the bottom of the resolved soil crust
				double WTDOLD = WTD;
				double SMCEQDEEP = SMCMAX * Math.Pow(-cell.PSISAT / (-cell.PSISAT - DZSNSO[NSOIL]), 1 / cell.BEXP);
				if (SMCWTD > SMCEQDEEP) {
					WTD = Math.Min((SMCWTD * DZSNSO[NSOIL]
					- SMCEQDEEP * ZSOIL0[NSOIL] + SMCMAX * (ZSOIL0[NSOIL] - DZSNSO[NSOIL])) /
					(SMCMAX - SMCEQDEEP), ZSOIL0[NSOIL]);
					RECH = -(WTDOLD - WTD) * (SMCMAX - SMCEQDEEP);
				} else {
					RECH = -(WTDOLD - (ZSOIL0[NSOIL] - DZSNSO[NSOIL])) * (SMCMAX - SMCEQDEEP);
					WTDOLD = ZSOIL0[NSOIL] - DZSNSO[NSOIL];
					//and now even further down
					double DZUP = (SMCEQDEEP - SMCWTD) * DZSNSO[NSOIL] / (SMCMAX - SMCEQDEEP);
					WTD = WTDOLD - DZUP;
					RECH -= (SMCMAX - SMCEQDEEP) * DZUP;
					SMCWTD = SMCEQDEEP;
				}
			}

			if (IWTD < NSOIL)
				SMCWTD = SMCMAX;
		}
		/// <summary>
		/// compute snow thermal conductivity and heat capacity
		/// 计算雪的导热(DF)与热容量HCPCT
		/// </summary>
		public static void ThermoProp(GridCell cell, int NSOIL, int NSNOW, int ISNOW, int IST, FortDoubleArray DZSNSO,
			double DT, double SNOWH, FortDoubleArray SNICE, FortDoubleArray SNLIQ, double CSOIL,
			FortDoubleArray SMC, FortDoubleArray SH2O, double TG, FortDoubleArray STC, double UR,
			double LAT, double Z0M, double ZLVL, int VEGTYP, int ISURBAN,
			FortDoubleArray DF, FortDoubleArray HCPCT, FortDoubleArray SNICEV, FortDoubleArray SNLIQV, FortDoubleArray EPORE,
			FortDoubleArray  FACT)
		{
			FortDoubleArray TKSNO = new FortDoubleArray(1 - Driver.NSnow, 0);
			FortDoubleArray CVSNO = new FortDoubleArray(1 - Driver.NSnow, 0);
			FortDoubleArray SICE = new FortDoubleArray(1, Driver.NSoil);
			CSNOW(ISNOW, Driver.NSnow, Driver.NSoil, SNICE, SNLIQ, DZSNSO,  //#input
				TKSNO, CVSNO, SNICEV, SNLIQV, EPORE);

			for (int IZ = ISNOW + 1; IZ <= 0; IZ++) {
				DF[IZ] = TKSNO[IZ];
				HCPCT[IZ] = CVSNO[IZ];
			}

			//compute soil thermal properties

			//这个索引没任何问题
			for (int IZ = 1; IZ <= Driver.NSoil; IZ++) {
				SICE[IZ] = SMC[IZ] - SH2O[IZ];
				HCPCT[IZ] = SH2O[IZ] * NoahMP.CWAT + (1.0 - cell.SMCMAX) * NoahMP.CSOIL + (cell.SMCMAX - SMC[IZ]) * NoahMP.CPAIR + SICE[IZ] * NoahMP.CICE;
				DF[IZ] = TDFCND(cell, SMC[IZ], SH2O[IZ]);
				if (DF[IZ] > 1000)
					throw new Exception();
			}

			if (VEGTYP == ISURBAN) {
				for (int IZ = 1; IZ <= Driver.NSoil; IZ++)
					DF[IZ] = 3.24;
			}
			// heat flux reduction effect from the overlying green canopy, adapted from
			// section 2.1.2 of Peters-Lidard et al. (1997, JGR, VOL 102(D4)).
			// not in use because of the separation of the canopy layer from the ground.
			// but this may represent the effects of leaf litter (Niu comments)
			//       DF1 = DF1 * Math.Exp (SBETA * SHDFAC)

			// compute lake thermal properties
			// (no consideration of turbulent mixing for this version)

			if (IST == 2) {
				//遇到Nsoil迭代变量的7元数组，需要使用IZ+Driver.NSnow,0,1,2层为雪层，3，4，5，6为土壤层
				for (int IZ = 1; IZ <= Driver.NSoil; IZ++) {
					if (STC[IZ] > NoahMP.TFRZ) {
						HCPCT[IZ] = NoahMP.CWAT;
						DF[IZ] = NoahMP.TKWAT;  //+ KEDDY * CWAT
					} else {
						HCPCT[IZ] = NoahMP.CICE;
						DF[IZ] = NoahMP.TKICE;
					}
					
				}
			}
			

			// combine a temporary variable used for melting/freezing of snow and frozen soil

			for (int IZ = ISNOW + 1; IZ <= Driver.NSoil; IZ++) {
				FACT[IZ] = DT / (HCPCT[IZ] * DZSNSO[IZ]);
//				if (Math.Abs(HCPCT[IZ]) < 1e-10 || Math.Abs(DZSNSO[IZ]) < 1e-10)
//					FACT[IZ] = 0;
				if (double.IsInfinity(DF[IZ]) || double.IsNaN(DF[IZ]))
					throw new Exception("");
			}
			// snow/soil interface
//			for (int IZ = 0; IZ < Driver.NSnow + Driver.NSoil; IZ++) {
//
//			if(double.IsNaN( DF[IZ]))
//				                    throw new Exception("");
//			}

			//如果雪线位于雪层顶部，说明几乎无雪
			if (ISNOW == 0) {
				DF[1] = (DF[1] * DZSNSO[1] + 0.35 * SNOWH) / (SNOWH + DZSNSO[1]);
			} else {
				DF[1] = (DF[1] * DZSNSO[1] + DF[0] * DZSNSO[0]) / (DZSNSO[0] + DZSNSO[1]);
			}
			for (int IZ = ISNOW + 1; IZ <= Driver.NSoil; IZ++) {
				
				if (double.IsNaN(DF[IZ]) || DF[IZ] > 1000)
					throw new Exception("");
			}
		}
		public static	void ESAT(double T, out double ESW, out double ESI, out double DESW, out double DESI)
		{
			//---------------------------------------------------------------------------------------------------
			// use polynomials to calculate saturation vapor pressure and derivative with
			// respect to temperature: over water when t > 0 c and over ice when t <= 0 c
			//---------------------------------------------------------------------------------------------------
			// in

			//  REAL, intent(in)  :: T              //temperature

			//out

			//  REAL, intent(out) :: ESW            //saturation vapor pressure over water (pa)
			//  REAL, intent(out) :: ESI            //saturation vapor pressure over ice (pa)
			//  REAL, intent(out) :: DESW           //d(esat)/dt over water (pa/K)
			//  REAL, intent(out) :: DESI           //d(esat)/dt over ice (pa/K)

			// local

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
		public static double TDC(double T)
		{
			return Math.Min(50, Math.Max(-50, T - NoahMP.TFRZ));
			//return Math.Min(30, Math.Max(-30, T - NoahMP.TFRZ));  //刘永和改成了上下30
		}
		/// <summary>
		/// 裸土通量
		/// </summary>
		static void BARE_FLUX(int NSNOW, int NSOIL, int ISNOW, double DT, double SAG, double //in
		               LWDN, double UR, double UU, double VV, double SFCTMP, double  //in
		               THAIR, double QAIR, double EAIR, double RHOAIR, double SNOWH, FortDoubleArray  //in
		               DZSNSO, double ZLVL, double ZPD, double Z0M, double           //in
		               EMG, FortDoubleArray STC, FortDoubleArray DF, double RSURF, double LATHEA,  //in
			double GAMMA, double RHSUR, int ILOC, int JLOC, double Q2,  //in
			ref double TGB, ref double CM, ref double CH,           //inout
			out double TAUXB, out double TAUYB, out double IRB, out double SHB, out double EVB,   //out
			out double  GHB, out double T2MB, double DX, double DZ8W, double IVGTYP, //out
			double QC, double PBLH, ref double QSFC, double PSFC, int ISURBAN,   //in
			int IZ0TLND, double SFCPRS, out double Q2B, out double EHB2)      //in
		{

			double EHB = 0;    //bare ground heat conductance
			double U10B = 0;    //10 m wind speed in eastward dir (m/s)
			double V10B = 0;    //10 m wind speed in eastward dir (m/s)
			double WSPD = 0;
			double WSTAR = 0;      //friction velocity n vertical direction (m/s) (only for SFCDif2)
			double Z0H = 0;        //roughness length, sensible heat, ground (m)
			double CSH = 0;        //coefficients for sh as function of ts
			double CEV = 0;        //coefficients for ev as function of esat[ts]
			double CH2B = 0;       //exchange coefficient for 2m temp.
			double CQ2B = 0;       //exchange coefficient for 2m temp.
			int VEGTYP = 0;     //vegetation type set to isbarren
			double ESTG = 0;       //saturation vapor pressure at tg (pa)
			double DESTG = 0;      //d(es)/dt at tg (pa/K)
			double ESATW = 0;      //es for water
			double ESATI = 0;      //es for ice
			double DSATW = 0;      //d(es)/dt at tg (pa/K) for water
			double DSATI = 0;      //d(es)/dt at tg (pa/K) for ice
			double MOZ = 0;        //Monin-Obukhov stability parameter
			double MOZOLD = 0;     //Monin-Obukhov stability parameter from prior iteration
			double FM = 0;         //momentum stability correction, weighted by prior iters
			double FH = 0;         //sen heat stability correction, weighted by prior iters
			double MOZSGN = 0;  //number of times MOZ changes sign
			double FM2 = 0;          //Monin-Obukhov momentum adjustment at 2m
			double FH2 = 0;          //Monin-Obukhov heat adjustment at 2m
			double CH2 = 0;          //Surface exchange at 2m
			int NITERB = 5;
			// -----------------------------------------------------------------
			// initialization variables that do not depend on stability iteration
			// -----------------------------------------------------------------
			double MPE = 1E-6;
//			double DTG = 0;
			//double MOZSGN = 0;
			//double MOZOLD = 0;
			double H = 0;
			double QFX = 0;
			double FV = 0.1;
			SHB = 0;
			EVB = 0;
			IRB = 0;
			GHB = 0;
			T2MB = 0;
			EHB2 = 0;
			Q2B = 0;
			double CIR = EMG * NoahMP.SB;
			
			double CGH = 2 * DF[ISNOW + 1] / DZSNSO[ISNOW + 1];

//			if (ISNOW == 0)
//				Console.WriteLine();
			for (int ITER = 1; ITER <= NITERB; ITER++) {
				// begin stability iteration
				if (ITER == 1)
					Z0H = Z0M;
				else
					Z0H = Z0M; //* Math.Exp(-CZIL*0.4*258.2*SQRT(FV*Z0M))
				

				if (NoahMP.OPT_SFC == 1) {
//					if (double.IsNaN(H))
//						throw new Exception("");
					
					SFCDIF1.SFCDIF(ITER, SFCTMP, RHOAIR, H, QAIR, //in
						ZLVL, ZPD, Z0M, Z0H, UR, //in
						MPE, ILOC, JLOC,                 //in
						ref MOZ, ref MOZSGN, ref FM, ref FH, ref FM2, ref FH2, //inout
						out CM, out CH, ref FV, out CH2);          //out
//					if(TGB<0)
//						throw new Exception("");
				}

				if (NoahMP.OPT_SFC == 2) {
					SFCDIF2.SFCDIF(ITER, Z0M, TGB, THAIR, UR, //in
						NoahMP.CZIL, ZLVL, ILOC, JLOC,         //in
						ref CM, ref CH, ref MOZ, ref WSTAR,         //inout
						ref FV);                                   //out
					// Undo the multiplication by windspeed that SFCDif2
					// applies to exchange coefficients CH and CM:
					CH /= UR;
					CM /= UR;
					if (SNOWH > 0) {
						CM = Math.Min(0.01, CM);   // CM  CH are too large, causing
						CH = Math.Min(0.01, CH);   // computational instability
					}

				}
				if (NoahMP.OPT_SFC == 3) {
					VEGTYP = Driver.vegparams.ISBARREN;
					//没有测试完
					SFCDIF3.SFCDIF(ILOC, JLOC, TGB, QSFC, PSFC,  //in
						PBLH, Z0M, Z0M, VEGTYP, ISURBAN,  //in
						IZ0TLND, UR, ITER, NITERB, SFCTMP,  //in
						THAIR, QAIR, QC, ZLVL,          //in
						SFCPRS, ref FV, ref CM, ref CH, out CH2B,  //inout
						out CQ2B, out MOZ);                               //out
					// Undo the multiplication by windspeed that SFCDIF3
					// applies to exchange coefficients CH and CM:
					//throw new Exception("not implemented");
					CH /= UR;
					CM /= UR;
					CH2B /= UR;

					if (SNOWH > 0) {     // jref: does this still count??
						CM = Math.Min(0.01, CM);   // CM  CH are too large, causing
						CH = Math.Min(0.01, CH);   // computational instability
						CH2B = Math.Min(0.01, CH2B);
						CQ2B = Math.Min(0.01, CQ2B);
					}
				}
				
				if (NoahMP.OPT_SFC == 4) {
//					if (double.IsNaN(TGB))
//						throw new Exception("");
					SFCDIF4.SFCDIF(ILOC, JLOC, UU, VV, SFCTMP,  //in
						SFCPRS, PSFC, PBLH, DX, Z0M,
						TGB, QAIR, ZLVL, IZ0TLND, ref QSFC,
						H, QFX, ref CM, ref CH, ref CH2B,
						ref CQ2B, ref MOZ, ref FV, out U10B, out V10B);
//					if (TGB < 0)
//						throw new Exception("");
					//Console.WriteLine("TGB=" + TGB);
//					if (double.IsNaN(TGB) || double.IsNaN(TG))
//						throw new Exception("");
					// Undo the multiplication by windspeed that SFCDIF4
					// applies to exchange coefficients CH and CM:
					CH /= UR;
					CM /= UR;
					CH2B /= UR;
					if (SNOWH > 0) {     // jref: does this still count??
						CM = Math.Min(0.01, CM);   // CM  CH are too large, causing
						CH = Math.Min(0.01, CH);   // computational instability
						CH2B = Math.Min(0.01, CH2B);
						CQ2B = Math.Min(0.01, CQ2B);
					}
				} //end opt_sfc 4

				double RAMB = Math.Max(1, 1 / (CM * UR));
				double RAHB = Math.Max(1, 1 / (CH * UR));
				;
				double RAWB = RAHB;

				//jref - variables for diagnostics
				double EMB = 1 / RAMB;
				EHB = 1 / RAHB;
				if (NoahMP.OPT_SFC == 3 || NoahMP.OPT_SFC == 4) {
					double RAHB2 = Math.Max(1, 1 / (CH2B * UR));
					EHB2 = 1 / RAHB2;
					CQ2B = EHB2;
				}

				double T = TDC(TGB);
				ESAT(T, out ESATW, out ESATI, out DSATW, out DSATI);
				if (T > 0) {
					ESTG = ESATW;
					DESTG = DSATW;
				} else {
					ESTG = ESATI;
					DESTG = DSATI;
				}

				CSH = RHOAIR * NoahMP.CPAIR / RAHB;
				CEV = RHOAIR * NoahMP.CPAIR / GAMMA / (RSURF + RAWB);

				// surface fluxes and dtg
				IRB = CIR * Math.Pow(TGB, 4) - EMG * LWDN;
				SHB = CSH * (TGB - SFCTMP);
				EVB = CEV * (ESTG * RHSUR - EAIR);
				GHB = CGH * (TGB - STC[ISNOW + 1]);

				double B = SAG - IRB - SHB - EVB - GHB;
				double A = 4 * CIR * Math.Pow(TGB, 3) + CSH + CEV * DESTG + CGH;
				double DTG = B / A;
				if (double.IsInfinity(DTG))
					throw new Exception();

				IRB += 4 * CIR * Math.Pow(TGB, 3) * DTG;
				SHB += CSH * DTG;
				EVB += CEV * DESTG * DTG;
				GHB += CGH * DTG;
				if (double.IsNaN(EVB))
					throw new Exception();
				// update ground surface temperature
				TGB += DTG;

				// for M-O length
				H = CSH * (TGB - SFCTMP);

				T = TDC(TGB);
				ESAT(T, out ESATW, out ESATI, out DSATW, out DSATI);
				if (T > 0)
					ESTG = ESATW;
				else
					ESTG = ESATI;
				QSFC = 0.622 * (ESTG * RHSUR) / (PSFC - 0.378 * (ESTG * RHSUR));
				QFX = (QSFC - QAIR) * CEV * GAMMA / NoahMP.CPAIR;
				
			} // end stability iteration
			// -----------------------------------------------------------------

			if (double.IsNaN(EVB))
				throw new Exception();
			// if snow on ground and TG > TFRZ: reset TG = TFRZ. reevaluate ground fluxes.
			if (NoahMP.OPT_STC == 1) {
				if (SNOWH > 0.05 && TGB > NoahMP.TFRZ) {
					TGB = NoahMP.TFRZ;
					IRB = CIR * Math.Pow(TGB, 4) - EMG * LWDN;
					SHB = CSH * (TGB - SFCTMP);
					EVB = CEV * (ESTG * RHSUR - EAIR);          //ESTG reevaluate ?
					GHB = SAG - (IRB + SHB + EVB);
				}
			}
			// wind stresses

			TAUXB = -RHOAIR * CM * UR * UU;
			TAUYB = -RHOAIR * CM * UR * VV;
			//jref:start; errors in original equation corrected.
			// 2m air temperature
			if (NoahMP.OPT_SFC == 1 || NoahMP.OPT_SFC == 2) {
				EHB2 = FV * NoahMP.VKC / Math.Log((2 + Z0H) / Z0H);
				EHB2 = FV * NoahMP.VKC / (Math.Log((2 + Z0H) / Z0H) - FH2);
				CQ2B = EHB2;
				if (EHB2 < 1E-5) {
					T2MB = TGB;
					Q2B = QSFC;
				} else {
					
					T2MB = TGB - SHB / (RHOAIR * NoahMP.CPAIR) * 1 / EHB2;
					Q2B = QSFC - EVB / (LATHEA * RHOAIR) * (1 / CQ2B + RSURF);
				}
				if (VEGTYP == Driver.vegparams.ISURBAN)
					Q2B = QSFC;
			}

			// myj consistent 2m temperature over bare soil
			if (NoahMP.OPT_SFC == 3 || NoahMP.OPT_SFC == 4)
			if (EHB2 < 1E-5) {
				T2MB = TGB;
				Q2B = QSFC;
			} else {
				T2MB = TGB - SHB / (RHOAIR * NoahMP.CPAIR * EHB2);
				Q2B = QSFC - QFX / (RHOAIR * CQ2B);
			}
			CH = EHB;			
		}
		/// <summary>
		/// 
		/// </summary>
		/// <param name="VEGTYP"></param>
		/// <param name="IST"></param>
		/// <param name="ISC"></param>
		/// <param name="ICE"></param>
		/// <param name="NSOIL"></param>
		/// <param name="SNEQVO"></param>
		/// <param name="SNEQV"></param>
		/// <param name="DT"></param>
		/// <param name="COSZ"></param>
		/// <param name="SNOWH"></param>
		/// <param name="TG"></param>
		/// <param name="TV"></param>
		/// <param name="FSNO"></param>
		/// <param name="QSNOW"></param>
		/// <param name="FWET"></param>
		/// <param name="ELAI"></param>
		/// <param name="ESAI"></param>
		/// <param name="SMC"></param>
		/// <param name="SOLAD"></param>
		/// <param name="SOLAI"></param>
		/// <param name="FVEG"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="ALBOLD"></param>
		/// <param name="TAUSS"></param>
		/// <param name="FSUN">sunlit fraction of canopy [-]</param>
		/// <param name="LAISUN">sunlit leaf area (-)</param>
		/// <param name="LAISHA">shaded leaf area (-)</param>
		/// <param name="PARSUN">average absorbed par for sunlit leaves (w/m2)</param>
		/// <param name="PARSHA">average absorbed par for shaded leaves (w/m2)</param>
		/// <param name="SAV"></param>
		/// <param name="SAG"></param>
		/// <param name="FSR"></param>
		/// <param name="FSA"></param>
		/// <param name="FSRV"></param>
		/// <param name="FSRG"></param>
		/// <param name="BGAP"></param>
		/// <param name="WGAP"></param>
		static void	Radiation(int VEGTYP, int IST, int  ISC, int ICE, int NSOIL,
			double SNEQVO, double SNEQV, double DT, double COSZ, double SNOWH,
			double TG, double TV, double FSNO, double QSNOW, double FWET,
			double ELAI, double ESAI, FortDoubleArray SMC, double[] SOLAD, double[] SOLAI,
			double FVEG, int ILOC, int JLOC,
			ref double ALBOLD, ref double TAUSS,                            // //inout
			out double FSUN, out double LAISUN, out double LAISHA, out double PARSUN, out double PARSHA, // //out
			out double SAV, out double SAG, out double FSR, out double FSA, out double FSRV, //
			out double  FSRG, out double BGAP, out double WGAP)
		{
			
			// surface abeldo
			double MPE = 1E-6;
			

			double[] ALBGRI, FTDD;
			double FSHA, VAI;
			double[] ALBD, ALBI, FABI, FABD, FTID, FTII;
			double[] ALBGRD = new double[2];
			ALBGRI = new double[2];
			ALBD = new double[2];
			ALBI = new double[2];
			FABD = new double[2];
			FABI = new double[2];
			FTDD = new double[2];
			FTID = new double[2];
			FTII = new double[2];
			double[] FREVD = new double[2];
			double[] FREVI = new double[2];
			double[] FREGD = new double[2];
			double[] FREGI = new double[2];
			int VEG = 0;
			
			
			//FREGD = 0;
			
			BGAP = 0;
			WGAP = 0;
			FSUN = 0; //sunlit fraction of canopy [-]
			FSRG = FSRV = 0;
			
			Albedo(VEGTYP, IST, ISC, ICE, Driver.NSoil,
				DT, COSZ, ELAI, ESAI,
				TG, TV, SNOWH, FSNO, FWET,
				SMC, SNEQVO, SNEQV, QSNOW, FVEG,
				ILOC, JLOC,
				ref ALBOLD, ref TAUSS, // //inout
				ALBGRD, ALBGRI, ALBD, ALBI, FABD, // //out
				FABI, FTDD, FTID, FTII, out FSUN, // //)   //out
				FREVI, FREVD, FREGD, FREGI, out BGAP, // //inout
				out WGAP);
			// surface radiation

			FSHA = 1 - FSUN;
			LAISUN = ELAI * FSUN;
			LAISHA = ELAI * FSHA;
			VAI = ELAI + ESAI;
			if (VAI > 0)
				VEG = 1;
			else
				VEG = 0;
			

			SURRAD(MPE, FSUN, FSHA, ELAI, VAI,
				LAISUN, LAISHA, SOLAD, SOLAI, FABD,
				FABI, FTDD, FTID, FTII, ALBGRD,
				ALBGRI, ALBD, ALBI, ILOC, JLOC,
				out PARSUN, out PARSHA, out SAV, out SAG, out FSA, // //out
				out FSR,                                //  //out
				FREVI, FREVD, FREGD, FREGI, ref FSRV, // //inout
				ref FSRG);
		}
		public static	void SURRAD(double MPE, double FSUN, double FSHA, double ELAI, double VAI,
			double LAISUN, double LAISHA, double[] SOLAD, double[] SOLAI, double[] FABD,
			double[] FABI, double[] FTDD, double[] FTID, double[] FTII, double[] ALBGRD,
			double[] ALBGRI, double[] ALBD, double[] ALBI, int ILOC, int JLOC,
			out double PARSUN, out double PARSHA, out double SAV, out double SAG, out double FSA, // //out
			out double FSR,                                //       //out
			double[]  FREVI, double[] FREVD, double[] FREGD, double[] FREGI, ref double FSRV,
			ref double  FSRG)
		{

			// --------------------------------------------------------------------------------------------------
			//  IMPLICIT NONE
			// --------------------------------------------------------------------------------------------------
			// input
//
			//  INTEGER, INTENT(IN)              :: ILOC
			//  INTEGER, INTENT(IN)              :: JLOC
			//  REAL, INTENT(IN)                 :: MPE     //prevents underflow errors if division by zero
//
			//  REAL, INTENT(IN)                 :: FSUN    //sunlit fraction of canopy
			//  REAL, INTENT(IN)                 :: FSHA    //shaded fraction of canopy
			//  REAL, INTENT(IN)                 :: ELAI    //leaf area, one-sided
			//  REAL, INTENT(IN)                 :: VAI     //leaf + stem area, one-sided
			//  REAL, INTENT(IN)                 :: LAISUN  //sunlit leaf area index, one-sided
			//  REAL, INTENT(IN)                 :: LAISHA  //shaded leaf area index, one-sided
//
			//  REAL, DIMENSION(1:2), INTENT(IN) :: SOLAD   //incoming direct solar radiation (w/m2)
			//  REAL, DIMENSION(1:2), INTENT(IN) :: SOLAI   //incoming diffuse solar radiation (w/m2)
			//  REAL, DIMENSION(1:2), INTENT(IN) :: FABD    //flux abs by veg (per unit incoming direct flux)
			//  REAL, DIMENSION(1:2), INTENT(IN) :: FABI    //flux abs by veg (per unit incoming diffuse flux)
			//REAL, DIMENSION(1:2), INTENT(IN) :: FTDD    //down dir flux below veg (per incoming dir flux)
			//  REAL, DIMENSION(1:2), INTENT(IN) :: FTID    //down dif flux below veg (per incoming dir flux)
			//  REAL, DIMENSION(1:2), INTENT(IN) :: FTII    //down dif flux below veg (per incoming dif flux)
			//  REAL, DIMENSION(1:2), INTENT(IN) :: ALBGRD  //ground albedo (direct)
			//  REAL, DIMENSION(1:2), INTENT(IN) :: ALBGRI  //ground albedo (diffuse)
			//  REAL, DIMENSION(1:2), INTENT(IN) :: ALBD    //overall surface albedo (direct)
			//  REAL, DIMENSION(1:2), INTENT(IN) :: ALBI    //overall surface albedo (diffuse)
//
			//  REAL, DIMENSION(1:2), INTENT(IN) :: FREVD    //overall surface albedo veg (direct)
			//  REAL, DIMENSION(1:2), INTENT(IN) :: FREVI    //overall surface albedo veg (diffuse)
			//  REAL, DIMENSION(1:2), INTENT(IN) :: FREGD    //overall surface albedo grd (direct)
			//  REAL, DIMENSION(1:2), INTENT(IN) :: FREGI    //overall surface albedo grd (diffuse)
//
			// output
//
			//  REAL, INTENT(OUT)                :: PARSUN  //average absorbed par for sunlit leaves (w/m2)
			//  REAL, INTENT(OUT)                :: PARSHA  //average absorbed par for shaded leaves (w/m2)
			//  REAL, INTENT(OUT)                :: SAV     //solar radiation absorbed by vegetation (w/m2)
			//  REAL, INTENT(OUT)                :: SAG     //solar radiation absorbed by ground (w/m2)
			//  REAL, INTENT(OUT)                :: FSA     //total absorbed solar radiation (w/m2)
			//  REAL, INTENT(OUT)                :: FSR     //total reflected solar radiation (w/m2)
			//  REAL, INTENT(OUT)                :: FSRV    //reflected solar radiation by vegetation
			//  REAL, INTENT(OUT)                :: FSRG    //reflected solar radiation by ground
//
//			// ------------------------ local variables ----------------------------------------------------
			//  INTEGER                          :: IB      //waveband number (1=vis, 2=nir)
			//  INTEGER                          :: NBAND   //number of solar radiation waveband classes
//
			//  REAL                             :: ABS     //absorbed solar radiation (w/m2)
			//  REAL                             :: RNIR    //reflected solar radiation [nir] (w/m2)
			//  REAL                             :: RVIS    //reflected solar radiation [vis] (w/m2)
			//  REAL                             :: LAIFRA  //leaf area fraction of canopy
			//  REAL                             :: TRD     //transmitted solar radiation: direct (w/m2)
			//  REAL                             :: TRI     //transmitted solar radiation: diffuse (w/m2)
			//  REAL, DIMENSION(1:2)             :: CAD     //direct beam absorbed by canopy (w/m2)
			//  REAL, DIMENSION(1:2)             :: CAI     //diffuse radiation absorbed by canopy (w/m2)
			// ---------------------------------------------------------------------------------------------

			int NBAND = 2;

			// zero summed solar fluxes

			SAG = 0;
			SAV = 0;
			FSA = 0;

			// loop over nband wavebands
			double[] CAD = new double[NBAND];
			double[] CAI = new double[NBAND];
			for (int IB = 0; IB < NBAND; IB++) {

				// absorbed by canopy
				CAD[IB] = SOLAD[IB] * FABD[IB];
				CAI[IB] = SOLAI[IB] * FABI[IB];
				SAV += CAD[IB] + CAI[IB];
				FSA += CAD[IB] + CAI[IB];

				// transmitted solar fluxes incident on ground

				double TRD = SOLAD[IB] * FTDD[IB];
				double TRI = SOLAD[IB] * FTID[IB] + SOLAI[IB] * FTII[IB];

				// solar radiation absorbed by ground surface

				double ABS = TRD * (1 - ALBGRD[IB]) + TRI * (1 - ALBGRI[IB]);
				if (double.IsNaN(ABS))
					throw new Exception();
				SAG += ABS;
				FSA += ABS;
			}

			// partition visible canopy absorption to sunlit and shaded fractions
			// to get average absorbed par for sunlit and shaded leaves

			double LAIFRA = ELAI / Math.Max(VAI, MPE);
			if (FSUN > 0) {
				PARSUN = (CAD[0] + FSUN * CAI[0]) * LAIFRA / Math.Max(LAISUN, MPE);
				PARSHA = (FSHA * CAI[0]) * LAIFRA / Math.Max(LAISHA, MPE);
				
			} else {
				PARSUN = 0;
				PARSHA = (CAD[0] + CAI[0]) * LAIFRA / Math.Max(LAISHA, MPE);
			}
//			if (double.IsNaN(PARSHA))
//				throw new Exception("");
			// reflected solar radiation

			double RVIS = ALBD[0] * SOLAD[0] + ALBI[0] * SOLAI[0];
			double RNIR = ALBD[1] * SOLAD[1] + ALBI[1] * SOLAI[1];
			FSR = RVIS + RNIR;

			// reflected solar radiation of veg. and ground (combined ground)
			FSRV = FREVD[0] * SOLAD[0] + FREVI[0] * SOLAI[0] + FREVD[1] * SOLAD[1] + FREVI[1] * SOLAI[1];
			FSRG = FREGD[0] * SOLAD[0] + FREGI[0] * SOLAI[0] + FREGD[1] * SOLAD[1] + FREGI[1] * SOLAI[1];

		}
		/// <summary>
		/// 已调过，应该没问题了
		/// </summary>
		/// <param name="ITER"></param>
		/// <param name="VAI"></param>
		/// <param name="RHOAIR">density air (kg/m3)</param>
		/// <param name="HG">sensible heat flux, 感热</param>
		/// <param name="TAH">air temperature at height z0h+zpd (k)</param>
		/// <param name="ZPD">zero plane displacement (m)</param>
		/// <param name="Z0MG"></param>
		/// <param name="Z0HG"></param>
		/// <param name="HCAN"></param>
		/// <param name="UC"></param>
		/// <param name="Z0H"></param>
		/// <param name="FV"></param>
		/// <param name="CWP"></param>
		/// <param name="VEGTYP"></param>
		/// <param name="MPE"></param>
		/// <param name="TV"></param>
		/// <param name="MOZG"></param>
		/// <param name="FHG"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="RAMG"></param>
		/// <param name="RAHG"></param>
		/// <param name="RAWG"></param>
		/// <param name="RB"></param>
		public static	void RAGRB(int ITER, double VAI, double RHOAIR, double HG, double TAH, double   //in
		                  ZPD, double Z0MG, double Z0HG, double HCAN, double UC, double   //in
		                  Z0H, double FV, double CWP, int VEGTYP, double MPE, double   //in
		                  TV, ref double MOZG, ref double FHG, int ILOC, int JLOC,  //inout
			out double  RAMG, out double RAHG, out double RAWG, out double RB)           //out
		{
			// --------------------------------------------------------------------------------------------------
			// compute under-canopy aerodynamic resistance RAG and leaf boundary layer
			// resistance RB
			// --------------------------------------------------------------------------------------------------
			//  USE NoahMP3_VEG_PARAMETERS
			// --------------------------------------------------------------------------------------------------
			// --------------------------------------------------------------------------------------------------
			// inputs

			//  INTEGER,              INTENT(IN) :: ILOC   //grid index
			//  INTEGER,              INTENT(IN) :: JLOC   //grid index
			//  INTEGER,              INTENT(IN) :: ITER   //iteration index
			//  INTEGER,              INTENT(IN) :: VEGTYP //vegetation physiology type
			//  REAL,                 INTENT(IN) :: VAI    //total LAI + stem area index, one sided
			//  REAL,                 INTENT(IN) :: RHOAIR //density air (kg/m3)
			//  REAL,                 INTENT(IN) :: HG     //ground sensible heat flux (w/m2)
			//  REAL,                 INTENT(IN) :: TV     //vegetation temperature (k)
			//  REAL,                 INTENT(IN) :: TAH    //air temperature at height z0h+zpd (k)
			//  REAL,                 INTENT(IN) :: ZPD    //zero plane displacement (m)
			//  REAL,                 INTENT(IN) :: Z0MG   //roughness length, momentum, ground (m)
			//  REAL,                 INTENT(IN) :: HCAN   //canopy height (m) [note: hcan >= z0mg]
			//  REAL,                 INTENT(IN) :: UC     //wind speed at top of canopy (m/s)
			//  REAL,                 INTENT(IN) :: Z0H    //roughness length, sensible heat (m)
			//  REAL,                 INTENT(IN) :: Z0HG   //roughness length, sensible heat, ground (m)
			//  REAL,                 INTENT(IN) :: FV     //friction velocity (m/s)
			//  REAL,                 INTENT(IN) :: CWP    //canopy wind parameter
			//  REAL,                 INTENT(IN) :: MPE    //prevents overflow error [ZZ] division by zero
			// in  out
//
			//  REAL,              INTENT(INOUT) :: MOZG   //Monin-Obukhov stability parameter
			//  REAL,              INTENT(INOUT) :: FHG    //stability correction
//
			// outputs
			//  REAL                             :: RAMG   //aerodynamic resistance for momentum (s/m)
			//  REAL                             :: RAHG   //aerodynamic resistance for sensible heat (s/m)
			//  REAL                             :: RAWG   //aerodynamic resistance for water vapor (s/m)
			//  REAL                             :: RB     //bulk leaf boundary layer resistance (s/m)
//
//
			//  REAL :: KH           //turbulent transfer coefficient, sensible heat, (m2/s)
			double TMP1 = 0;        //temporary calculation
			double TMP2 = 0;         //temporary calculation
			//  REAL :: TMPRAH2      //temporary calculation for aerodynamic resistances
			//  REAL :: TMPRB        //temporary calculation for rb
			//  real :: MOLG,FHGNEW,CWPC
			// --------------------------------------------------------------------------------------------------
			// stability correction to below canopy resistance

			MOZG = 0;
			double MOLG = 0;

			if (ITER > 1) {
				TMP1 = NoahMP.VKC * (NoahMP.GRAV / TAH) * HG / (RHOAIR * NoahMP.CPAIR);
				if (Math.Abs(TMP1) <= MPE)
					TMP1 = MPE;
				MOLG = -1 * FV * FV * FV / TMP1;
				MOZG = Math.Min((ZPD - Z0MG) / MOLG, 1);
			}

			double FHGNEW = 0;
			if (MOZG < 0) {
				FHGNEW = Math.Pow(1 - 15 * MOZG, -0.25);
			} else {
				FHGNEW = 1 + 4.7 * MOZG;
			}

			if (ITER == 1)
				FHG = FHGNEW;
			else
				FHG = 0.5 * (FHG + FHGNEW);
			

			if (HCAN < 0.00001)
				HCAN = 0.001;
			double CWPC = Math.Sqrt(CWP * VAI * HCAN * FHG);
			//       CWPC = (CWP*FHG)**0.5

			TMP1 = Math.Exp(-CWPC * Z0HG / HCAN);
			TMP2 = Math.Exp(-CWPC * (Z0H + ZPD) / HCAN);
			double TMPRAH2 = HCAN * Math.Exp(CWPC) / CWPC * (TMP1 - TMP2);

			// aerodynamic resistances raw and rah between heights zpd+z0h and z0hg.
			double KH = Math.Max(NoahMP.VKC * FV * (HCAN - ZPD), MPE);
			RAMG = 0;
			RAHG = TMPRAH2 / KH;
			RAWG = RAHG;

			// leaf boundary layer resistance

			double TMPRB = CWPC * 50 / (1 - Math.Exp(-CWPC / 2));
			RB = TMPRB * Math.Sqrt(Driver.vegparams.DLEAF[VEGTYP - 1] / UC);
			//       RB = 200
		}
		/// <summary>
		/// 2022年11月16日又检查了多遍
		/// </summary>
		/// <param name="NSOIL"></param>
		/// <param name="ZSOIL"></param>
		/// <param name="DT"></param>
		/// <param name="PDDUM"></param>
		/// <param name="ETRANI"></param>
		/// <param name="QSEVA"></param>
		/// <param name="SH2O"></param>
		/// <param name="SMC"></param>
		/// <param name="ZWT"></param>
		/// <param name="FCR"></param>
		/// <param name="SICEMAX"></param>
		/// <param name="FCRMAX"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="SMCWTD"></param>
		/// <param name="RHSTT"></param>
		/// <param name="AI"></param>
		/// <param name="BI"></param>
		/// <param name="CI"></param>
		/// <param name="QDRAIN"></param>
		/// <param name="WCND">hydraulic conductivity, AI,BI,CI的初次生成是在这个函数</param>
		public static void SRT(GridCell cell, int NSOIL, FortDoubleArray ZSOIL, double DT, double PDDUM, FortDoubleArray ETRANI,
			double QSEVA, FortDoubleArray SH2O, FortDoubleArray SMC, double ZWT, FortDoubleArray FCR, //in
			double  SICEMAX, double FCRMAX, int ILOC, int JLOC, double SMCWTD,  //in
			FortDoubleArray RHSTT, FortDoubleArray AI, FortDoubleArray BI, FortDoubleArray CI, out double QDRAIN, FortDoubleArray WCND)  //out
		{
			// ----------------------------------------------------------------------
			// 

//			input

			//    INTEGER,                  INTENT(IN)   ILOC   //grid index
			//    INTEGER,                  INTENT(IN)   JLOC   //grid index
			//    INTEGER,                  INTENT(IN)   NSOIL
			//    double[], INTENT(IN)   ZSOIL
			//    REAL,                     INTENT(IN)   DT
			//    REAL,                     INTENT(IN)   PDDUM
			//    REAL,                     INTENT(IN)   QSEVA
			//    double[], INTENT(IN)   ETRANI
			//    double[], INTENT(IN)   SH2O
			//    double[], INTENT(IN)   SMC
			//    REAL,                     INTENT(IN)   ZWT    // water table depth [m]
			//    double[], INTENT(IN)   FCR
			//    REAL, INTENT(IN)                       FCRMAX //maximum of FCR (-)
			//    REAL,                     INTENT(IN)   SICEMAX//maximum soil ice content (m3/m3)
			//    REAL,                     INTENT(IN)   SMCWTD //soil moisture between bottom of the soil and the water table
//
			//// output
//
			//    double[], INTENT(OUT)  RHSTT
			//    double[], INTENT(OUT)  AI
			//    double[], INTENT(OUT)  BI
			//    double[], INTENT(OUT)  CI
			//    double[], INTENT(OUT)  WCND    //
			//    REAL,     INTENT(OUT)  QDRAIN  //
//			 
			// local
			int K;
			FortDoubleArray DDZ = new FortDoubleArray(1, Driver.NSoil); //FortDoubleArray
			FortDoubleArray DENOM = new FortDoubleArray(1, Driver.NSoil);
			FortDoubleArray DSMDZ = new FortDoubleArray(1, Driver.NSoil);
			FortDoubleArray WFLUX = new FortDoubleArray(1, Driver.NSoil);
			FortDoubleArray WDF = new FortDoubleArray(1, Driver.NSoil);
			FortDoubleArray SMX = new FortDoubleArray(1, Driver.NSoil);
			double TEMP1 = 0;
			double SMXWTD = 0; //soil moisture between bottom of the soil and water table
			double SMXBOT = 0;  //soil moisture below bottom to calculate flux

			// Niu and Yang (2006), J. of Hydrometeorology
			// ----------------------------------------------------------------------
			QDRAIN = 0;
			if (NoahMP.OPT_INF == 1) {  //最开始测试的就是该选项，没问题
				for (K = 1; K <= NSOIL; K++) {
					double wdf = WDF[K];
					double wcnd = WCND[K];
					
					WDFCND1(cell, out wdf, out wcnd, SMC[K], FCR[K]);
					WDF[K] = wdf;
					WCND[K] = wcnd;
					SMX[K] = SMC[K];
				}
				if (NoahMP.OPT_RUN == 5)
					SMXWTD = SMCWTD; //soil moisture between bottom of the soil and water table
			}
			if (NoahMP.OPT_INF == 2) {
				for (K = 1; K <= NSOIL; K++) {
					double wdf = WDF[K];
					double wcnd = WCND[K];
					WDFCND2(cell, out wdf, out wcnd, SH2O[K], SICEMAX); //调试时修正了一个错误
					WDF[K] = wdf;
					WCND[K] = wcnd;
					SMX[K] = SH2O[K];
				}

				if (NoahMP.OPT_RUN == 5)
					SMXWTD = SMCWTD * SH2O[NSOIL] / SMC[NSOIL];  //same liquid fraction as in the bottom layer
			}

//			Console.WriteLine("WCND=" + WCND[4]);
			for (K = 1; K <= NSOIL; K++) {
				if (K == 1) {
					DENOM[K] = -ZSOIL[K]; //当前土壤深度
					TEMP1 = -ZSOIL[K + 1];  //下一层土壤深度
					DDZ[K] = 2.0 / TEMP1;
					DSMDZ[K] = 2.0 * (SMX[K] - SMX[K + 1]) / TEMP1;
					WFLUX[K] = WDF[K] * DSMDZ[K] + WCND[K] - PDDUM + ETRANI[K] + QSEVA;
				} else if (K < NSOIL) {
					DENOM[K] = (ZSOIL[K - 1] - ZSOIL[K]);
					TEMP1 = (ZSOIL[K - 1] - ZSOIL[K + 1]);
					DDZ[K] = 2.0 / TEMP1;
					DSMDZ[K] = 2.0 * (SMX[K] - SMX[K + 1]) / TEMP1;
					WFLUX[K] = WDF[K] * DSMDZ[K] + WCND[K] - WDF[K - 1] * DSMDZ[K - 1] - WCND[K - 1] + ETRANI[K];
				} else {  //K==NSOIL
					DENOM[K] = ZSOIL[K - 1] - ZSOIL[K]; //DENOM为第K层土壤的深度
					if (NoahMP.OPT_RUN == 1 || NoahMP.OPT_RUN == 2) {
						QDRAIN = 0;
					}
					if (NoahMP.OPT_RUN == 3) {
						//free drainage, SLOPE参数的影响比较大
						QDRAIN = cell.SLOPE * WCND[K];
//						if (QDRAIN > 1e-7)
//							throw new Exception();
					}
					if (NoahMP.OPT_RUN == 4) {
						QDRAIN = (1.0 - FCRMAX) * WCND[K];
					}
					if (NoahMP.OPT_RUN == 5) {   //gmm new m-m&f water table dynamics formulation
						TEMP1 = 2.0 * DENOM[K];
						if (ZWT < ZSOIL[NSOIL] - DENOM[NSOIL]) {
							//gmm interpolate from below, midway to the water table, to the middle of the auxiliary layer below the soil bottom
							SMXBOT = SMX[K] - (SMX[K] - SMXWTD) * DENOM[K] * 2 / (DENOM[K] + ZSOIL[K] - ZWT);
						} else {
							SMXBOT = SMXWTD;
						}
						DSMDZ[K] = 2.0 * (SMX[K] - SMXBOT) / TEMP1;
						QDRAIN = WDF[K] * DSMDZ[K] + WCND[K];
					}
					WFLUX[K] = -(WDF[K - 1] * DSMDZ[K - 1]) - WCND[K - 1] + ETRANI[K] + QDRAIN;
				}
			}  //end for

			for (K = 1; K <= NSOIL; K++) {
				if (K == 1) {
					AI[K] = 0.0;
					BI[K] = WDF[K] * DDZ[K] / DENOM[K];
					CI[K] = -BI[K];
				} else if (K < NSOIL) {
					AI[K] = -WDF[K - 1] * DDZ[K - 1] / DENOM[K];
					CI[K] = -WDF[K] * DDZ[K] / DENOM[K];
					BI[K] = -(AI[K] + CI[K]);
				} else {
					AI[K] = -WDF[K - 1] * DDZ[K - 1] / DENOM[K];
					CI[K] = 0.0;
					BI[K] = -(AI[K] + CI[K]);
				}

				RHSTT[K] = WFLUX[K] / (-DENOM[K]);				
			}
		}
		
		/// <summary>
		/// 计算气孔阻力和光合作用  canopy stomatal resistance
		/// </summary>
		/// <param name="VEGTYP"></param>
		/// <param name="MPE"></param>
		/// <param name="APAR"></param>
		/// <param name="FOLN"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="TV">foliage temperature (k)</param>
		/// <param name="EI">vapor pressure inside leaf (sat vapor press at tv) (pa)</param>
		/// <param name="EA">vapor pressure of canopy air (pa)</param>
		/// <param name="SFCTMP"></param>
		/// <param name="SFCPRS"></param>
		/// <param name="O2"></param>
		/// <param name="CO2"></param>
		/// <param name="IGS"></param>
		/// <param name="BTRAN">soil water transpiration factor (0 to 1)</param>
		/// <param name="RB">boundary layer resistance (s/m)</param>
		/// <param name="RS">leaf stomatal resistance (s/m)</param>
		/// <param name="PSN">光合作用量</param>
		public static void STOMATA(int VEGTYP, double MPE, double APAR, double FOLN, int ILOC, int JLOC,  //in
			double TV, double EI, double EA, double SFCTMP, double SFCPRS, double  //in
		                    O2, double CO2, double IGS, double BTRAN, double RB, //in
			out double RS, out double PSN)                              //out
		{
			// --------------------------------------------------------------------------------------------------
			//  USE NoahMP3_VEG_PARAMETERS

			// input
			//      INTEGER,INTENT(IN)  :: ILOC   //grid index
			//      INTEGER,INTENT(IN)  :: JLOC   //grid index
			//      INTEGER,INTENT(IN)  :: VEGTYP //vegetation physiology type
//
			//      REAL, INTENT(IN)    :: IGS    //growing season index (0=off, 1=on)
			//      REAL, INTENT(IN)    :: MPE    //prevents division by zero errors
//
			//      REAL, INTENT(IN)    :: TV     //foliage temperature (k)
			//      REAL, INTENT(IN)    :: EI     //vapor pressure inside leaf (sat vapor press at tv) (pa)
			//      REAL, INTENT(IN)    :: EA     //vapor pressure of canopy air (pa)
			//      REAL, INTENT(IN)    :: APAR   //par absorbed per unit lai (w/m2)
			//      REAL, INTENT(IN)    :: O2     //atmospheric o2 concentration (pa)
			//      REAL, INTENT(IN)    :: CO2    //atmospheric co2 concentration (pa)
			//      REAL, INTENT(IN)    :: SFCPRS //air pressure at reference height (pa)
			//      REAL, INTENT(IN)    :: SFCTMP //air temperature at reference height (k)
			//      REAL, INTENT(IN)    :: BTRAN  //soil water transpiration factor (0 to 1)
			//      REAL, INTENT(IN)    :: FOLN   //foliage nitrogen concentration (%)
			//      REAL, INTENT(IN)    :: RB     //boundary layer resistance (s/m)
//
			// output
			//      REAL, INTENT(OUT)   :: RS     //leaf stomatal resistance (s/m)
			//      REAL, INTENT(OUT)   :: PSN    //foliage photosynthesis (umol co2 /m2/ s) [always +]
			// in&out
			//      REAL                :: RLB    //boundary layer resistance (s m2 / umol)
			// ---------------------------------------------------------------------------------------------
//
			// ------------------------ local variables ----------------------------------------------------
			//      INTEGER :: ITER     //iteration index
			//      INTEGER :: NITER    //number of iterations
//
			int NITER = 3;
			//      SAVE NITER
//
			//      REAL :: AB          //used in statement functions
			//      REAL :: BC          //used in statement functions
			//      REAL :: F1          //generic temperature response (statement function)
			//      REAL :: F2          //generic temperature inhibition (statement function)
			//      REAL :: TC          //foliage temperature (degree Celsius)
			//      REAL :: CS          //co2 concentration at leaf surface (pa)
			//      REAL :: KC          //co2 Michaelis-Menten constant (pa)
			//      REAL :: KO          //o2 Michaelis-Menten constant (pa)
			//      REAL :: A,B,C,Q     //intermediate calculations for RS
			//      REAL :: R1,R2       //roots for RS
			//      REAL :: FNF         //foliage nitrogen adjustment factor (0 to 1)
			//      REAL :: PPF         //absorb photosynthetic photon flux (umol photons/m2/s)
			//      REAL :: WC          //Rubisco limited photosynthesis (umol co2/m2/s)
			//      REAL :: WJ          //light limited photosynthesis (umol co2/m2/s)
			//      REAL :: WE          //export limited photosynthesis (umol co2/m2/s)
			//      REAL :: CP          //co2 compensation point (pa)
			//      REAL :: CI          //internal co2 (pa)
			//      REAL :: AWC         //intermediate calculation for wc
			//      REAL :: VCMX        //maximum rate of carbonylation (umol co2/m2/s)
			//      REAL :: J           //electron transport (umol co2/m2/s)
//			 REAL :: CEA         //constrain ea or else model blows up
			//      REAL :: CF          //s m2/umol -> s/m


			//      REAL :: T
			// ---------------------------------------------------------------------------------------------

			// initialize RS=RSMAX and PSN=0 because will only do calculations
			// for APAR > 0, in which case RS <= RSMAX and PSN >= 0

			double CF = SFCPRS / (8.314 * SFCTMP) * 1e06;
			double C3PSN = Driver.vegparams.C3PSN[VEGTYP - 1];
			double MP = Driver.vegparams.MP[VEGTYP - 1];
			double BP = Driver.vegparams.BP[VEGTYP - 1];
			RS = 1 / BP * CF;
			PSN = 0;

			if (APAR <= 0)
				return;

			double FNF = Math.Min(FOLN / Math.Max(MPE, Driver.vegparams.FOLNMX[VEGTYP - 1]), 1.0); //FNF:foliage nitrogen adjustment factor (0 to 1)
			double TC = TV - NoahMP.TFRZ;
			double PPF = 4.6 * APAR;
			double J = PPF * Driver.vegparams.QE25[VEGTYP - 1];
			double KC = Driver.vegparams.KC25[VEGTYP - 1] * F1_Func(Driver.vegparams.AKC[VEGTYP - 1], TC);
			double KO = Driver.vegparams.KO25[VEGTYP - 1] * F1_Func(Driver.vegparams.AKO[VEGTYP - 1], TC);
			double AWC = KC * (1 + O2 / KO);
			double CP = 0.5 * KC / KO * O2 * 0.21;
			//maximum rate of carbonylation (umol co2/m2/s)
			double VCMX = Driver.vegparams.VCMX25[VEGTYP - 1] / F2_Func(TC) * FNF * BTRAN * F1_Func(Driver.vegparams.AVCMX[VEGTYP - 1], TC);

			// first guess ci			
			double CI = 0.7 * CO2 * C3PSN + 0.4 * CO2 * (1.0 - C3PSN); //CI: internal co2 (pa)
			// rb: s/m -> s m**2 / umol

			double	RLB = RB / CF;

			// constrain ea
			double CEA = Math.Max(0.25 * EI * C3PSN + 0.40 * EI * (1.0 - C3PSN), Math.Min(EA, EI));

			// ci iteration
			//jref: C3PSN is equal to 1 for all veg types.
			for (int ITER = 1; ITER <= NITER; ITER++) {
				double	WJ = Math.Max(CI - CP, 0) * J / (CI + 2 * CP) * C3PSN + J * (1.0 - C3PSN);
				double	WC = Math.Max(CI - CP, 0) * VCMX / (CI + AWC) * C3PSN + VCMX * (1.0 - C3PSN);
				double	WE = 0.5 * VCMX * C3PSN + 4000 * VCMX * CI / SFCPRS * (1.0 - C3PSN);
				PSN = Math.Min(Math.Min(WJ, WC), WE) * IGS; 
				
				double	CS = Math.Max(CO2 - 1.37 * RLB * SFCPRS * PSN, MPE);
				double	A = MP * PSN * SFCPRS * CEA / (CS * EI) + BP;
				double	B = (MP * PSN * SFCPRS / CS + BP) * RLB - 1;
				double	C = -RLB;
				double Q = -0.5 * (B - Math.Sqrt(B * B - 4.0 * A * C));
				if (B >= 0)
					Q = -0.5 * (B + Math.Sqrt(B * B - 4.0 * A * C));
				
				
				double R1 = Q / A;
				double R2 = C / Q;
				RS = Math.Max(R1, R2);  //leaf stomatal resistance (s/m)
				CI = Math.Max(CS - PSN * SFCPRS * 1.65 * RS, 0); //CI: internal co2 (pa)
			}
			// rs, rb:  s m**2 / umol -> s/m
//			Console.WriteLine("PSN="+PSN);
			RS *= CF;
		}
		/// <summary>
		/// 已调试过，应该没问题了
		/// </summary>
		/// <param name="NSNOW"></param>
		/// <param name="NSOIL"></param>
		/// <param name="ISNOW"></param>
		/// <param name="ZSNSO"></param>
		/// <param name="STC"></param>
		/// <param name="TBOT"></param>
		/// <param name="ZBOT"></param>
		/// <param name="DT"></param>
		/// <param name="DF"></param>
		/// <param name="HCPCT"></param>
		/// <param name="SSOIL"></param>
		/// <param name="PHI"></param>
		/// <param name="AI"></param>
		/// <param name="BI"></param>
		/// <param name="CI"></param>
		/// <param name="RHSTS"></param>
		/// <param name="BOTFLX"></param>
		public	static	 void HRT(int NSNOW, int NSOIL, int ISNOW, FortDoubleArray ZSNSO,
			FortDoubleArray  STC, double TBOT, double ZBOT, double DT,
			FortDoubleArray DF, FortDoubleArray HCPCT, double SSOIL, FortDoubleArray PHI,
			FortDoubleArray AI, FortDoubleArray BI, FortDoubleArray CI, FortDoubleArray RHSTS, out double  BOTFLX)
		{
			// ----------------------------------------------------------------------
			// ----------------------------------------------------------------------
			// calculate the right hand side of the time tendency term of the soil
			// thermal diffusion equation.  also to compute ( prepare ) the matrix
			// coefficients for the tri-diagonal matrix of the implicit time scheme.
			// input

			//    INTEGER,                         INTENT(IN)  :: NSOIL  //no of soil layers (4)
			//    INTEGER,                         INTENT(IN)  :: NSNOW  //maximum no of snow layers (3)
			//    INTEGER,                         INTENT(IN)  :: ISNOW  //actual no of snow layers
			//    REAL,                            INTENT(IN)  :: TBOT   //bottom soil temp. at ZBOT (k)
			//    REAL,                            INTENT(IN)  :: ZBOT   //depth of lower boundary condition (m)
			//                                                           //from soil surface not snow surface
			//    REAL,                            INTENT(IN)  :: DT     //time step (s)
			//    REAL,                            INTENT(IN)  :: SSOIL  //ground heat flux (w/m2)
			//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: ZSNSO  //depth of layer-bottom of snow/soil (m)
			//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: STC    //snow/soil temperature (k)
			//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: DF     //thermal conductivity [w/m/k]
			//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: HCPCT  //heat capacity [j/m3/k]
			//    REAL, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: PHI    //light through water (w/m2)

			// output
			//double[] RHSTS=new double[NSNOW+NSOIL];//right-hand side of the matrix
			//double[] AI=new double[NSNOW+NSOIL];  //left-hand side coefficient
			//double[] BI=new double[NSNOW+NSOIL];    //left-hand side coefficient
			//double[] CI=new double[NSNOW+NSOIL];     //left-hand side coefficient
			BOTFLX = 0; //energy influx from soil bottom (w/m2)  为最后一个输出变量

			// local

			//    INTEGER                                      :: K
			FortDoubleArray DDZ = new FortDoubleArray(1 - NSNOW, NSOIL);
			//    REAL, DIMENSION(-NSNOW+1:NSOIL)              :: DZ
			FortDoubleArray DENOM = new FortDoubleArray(1 - NSNOW, NSOIL);
			FortDoubleArray DTSDZ = new FortDoubleArray(1 - NSNOW, NSOIL);
			FortDoubleArray EFLUX = new FortDoubleArray(1 - NSNOW, NSOIL);
			//    REAL                                         :: TEMP1
			// ----------------------------------------------------------------------
			for (int K = ISNOW + 1; K <= NSOIL; K++) {
				double TEMP1 = 0;
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
			//上段代码完全没有问题
			for (int K = ISNOW + 1; K <= NSOIL; K++) {   //2022年11月16日又修正了一个错误
				if (K == ISNOW + 1) {
					AI[K] = 0.0;
					CI[K] = -DF[K] * DDZ[K] / DENOM[K];
					if (NoahMP.OPT_STC == 1) {
						BI[K] = -CI[K];
					}
					if (NoahMP.OPT_STC == 2) {
						BI[K] = -CI[K] + DF[K] / (0.5 * ZSNSO[K] * ZSNSO[K] * HCPCT[K]);
					}
				} else if (K < NSOIL) {  //2022年11月16日又修正了一个错误
					AI[K] = -DF[K - 1] * DDZ[K - 1] / DENOM[K];
					CI[K] = -DF[K] * DDZ[K] / DENOM[K];
					BI[K] = -(AI[K] + CI[K]);
				} else if (K == NSOIL) {  //2022年11月16日又修正了一个错误
					AI[K] = -DF[K - 1] * DDZ[K - 1] / DENOM[K];
					CI[K] = 0.0;
					BI[K] = -(AI[K] + CI[K]);
				}
				
				RHSTS[K] = EFLUX[K] / (-DENOM[K]);
			}
			
		}

		public static void SNOWFALL(int NSOIL, int NSNOW, double DT, double QSNOW, double SNOWHIN,
			double SFCTMP, int ILOC, int JLOC,
			ref int ISNOW, ref double SNOWH, FortDoubleArray DZSNSO, FortDoubleArray STC, FortDoubleArray SNICE,
			FortDoubleArray SNLIQ, ref double SNEQV)
		{
			// ----------------------------------------------------------------------
			// snow depth and density to account for the new snowfall.
			// new values of snow depth & density returned.

			// input
			/*
  int,                            INTENT(IN)  ILOC   //grid index
  int,                            INTENT(IN)  JLOC   //grid index
  int,                            INTENT(IN)  NSOIL  //no. of soil layers
  int,                            INTENT(IN)  NSNOW  //maximum no. of snow layers
  double,                               INTENT(IN)  DT     //main time step (s)
  double,                               INTENT(IN)  QSNOW  //snow at ground srf (mm/s) [+]
  double,                               INTENT(IN)  SNOWHIN//snow depth increasing rate (m/s)
  double,                               INTENT(IN)  SFCTMP //surface air temperature [k]

// input and output

  int,                         INTENT(INOUT)  ISNOW  //actual no. of snow layers
  double,                            INTENT(INOUT)  SNOWH  //snow depth [m]
  double,                            INTENT(INOUT)  SNEQV  //swow water equivalent [m]
double, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT)  DZSNSO //thickness of snow/soil layers (m)
  double, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT)  STC    //snow layer temperature [k]
  double, DIMENSION(-NSNOW+1:    0), INTENT(INOUT)  SNICE  //snow layer ice [mm]
  double, DIMENSION(-NSNOW+1:    0), INTENT(INOUT)  SNLIQ  //snow layer liquid water [mm]
			 */
			// local

			int NEWNODE;            // 0-no new layers, 1-creating new layers
			// ----------------------------------------------------------------------
			NEWNODE = 0;

			// shallow snow / no layer

			if (ISNOW == 0 && QSNOW > 0) {
				SNOWH = SNOWH + SNOWHIN * DT;
				SNEQV = SNEQV + QSNOW * DT;
			}

			// creating a new layer

			if (ISNOW == 0 && QSNOW > 0 && SNOWH >= 0.025) { //MB: change limit
				//    if(ISNOW == 0  && QSNOW>0. && SNOWH >= 0.05)
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
				SNICE[ISNOW + 1] += QSNOW * DT;
				DZSNSO[ISNOW + 1] += SNOWHIN * DT;
			}
		}
		/// 实际只更新DZSNSO
		/// </summary>
		/// <param name="NSNOW"></param>
		/// <param name="NSOIL"></param>
		/// <param name="DT"></param>
		/// <param name="STC"></param>
		/// <param name="SNICE"></param>
		/// <param name="SNLIQ"></param>
		/// <param name="ZSOIL"></param>
		/// <param name="IMELT"></param>
		/// <param name="FICEOLD"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="ISNOW">不更新</param>
		/// <param name="DZSNSO"></param>
		/// <param name="ZSNSO">不更新</param>
		public static void COMPACT(int NSNOW, int NSOIL, double DT, FortDoubleArray STC, FortDoubleArray SNICE,
			FortDoubleArray  SNLIQ, FortDoubleArray ZSOIL, FortIntArray IMELT, FortDoubleArray FICEOLD, int ILOC, int JLOC,
			int ISNOW, FortDoubleArray DZSNSO, FortDoubleArray ZSNSO)
		{

			// input
			/*
   int,                         INTENT(IN)     ILOC   //grid index
   int,                         INTENT(IN)     JLOC   //grid index
   int,                         INTENT(IN)     NSOIL  //no. of soil layers [ =4]
   int,                         INTENT(IN)     NSNOW  //maximum no. of snow layers [ =3]
   int, DIMENSION(-NSNOW+1:0) , INTENT(IN)     IMELT  //melting state index [0-no melt;1-melt]
   double,                            INTENT(IN)     DT     //time step (sec)
   double, DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)     STC    //snow layer temperature [k]
   double, DIMENSION(-NSNOW+1:    0), INTENT(IN)     SNICE  //snow layer ice [mm]
   double, DIMENSION(-NSNOW+1:    0), INTENT(IN)     SNLIQ  //snow layer liquid water [mm]
   double, DIMENSION(       1:NSOIL), INTENT(IN)     ZSOIL  //depth of layer-bottom from soil srf
   double, DIMENSION(-NSNOW+1:    0), INTENT(IN)     FICEOLD//ice fraction at last timestep

// input and output
   int,                         INTENT(IN)  ISNOW  // actual no. of snow layers
   double, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT)  DZSNSO // snow layer thickness [m]
   double, DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT)  ZSNSO  // depth of snow/soil layer-bottom
			 */
			// local
			double C2 = 21e-3; //[m3/kg] ; // default 21.e-3
			double C3 = 2.5e-6; //[1/s]
			double C4 = 0.04; //[1/k]
			double C5 = 2.0; //
			double DM = 100.0; //upper Limit on destructive metamorphism compaction [kg/m3]
			double ETA0 = 0.8e+6; //viscosity coefficient [kg-s/m2]
			; //according to Anderson, it is between 0.52e6~1.38e6
			double BURDEN; //pressure of overlying snow [kg/m2]
			double DDZ1; //rate of settling of snow pack due to destructive metamorphism.
			double DDZ2; //rate of compaction of snow pack due to overburden.
			double DDZ3; //rate of compaction of snow pack due to melt [1/s]
			double DEXPF; //EXPF=exp(-c4*(273.15-STC)).
			double TD; //STC - TFRZ [K]
			double PDZDTC; //nodal rate of change in fractional-thickness due to compaction [fraction/s]
			double VOID; //void (1 - SNICE - SNLIQ)
			double WX; //water mass (ice + liquid) [kg/m2]
			double BI; //partial density of ice [kg/m3]
			FortDoubleArray FICE = new FortDoubleArray(1 - NSNOW, 0); //fraction of ice at current time step

			
			// ----------------------------------------------------------------------
			BURDEN = 0.0;

			for (int J = ISNOW + 1; J <= 0; J++) {
				//DO J = ISNOW+1, 0

				WX = SNICE[J] + SNLIQ[J];
				FICE[J] = SNICE[J] / WX;
				VOID = 1 - (SNICE[J] / NoahMP.DENICE + SNLIQ[J] / NoahMP.DENH2O) / DZSNSO[J];

				// Allow compaction only for non-saturated node and higher ice lens node.
				if (VOID > 0.001 && SNICE[J] > 0.1) {

					BI = SNICE[J] / DZSNSO[J];
					TD = Math.Max(0, NoahMP.TFRZ - STC[J]);
					DEXPF = Math.Exp(-C4 * TD);

					// Settling as a result of destructive metamorphism

					DDZ1 = -C3 * DEXPF;

					if (BI > DM) {
						DDZ1 = DDZ1 * Math.Exp(-46.0E-3 * (BI - DM));
					}
					// Liquid water term

					if (SNLIQ[J] > 0.01 * DZSNSO[J]) {
						DDZ1 = DDZ1 * C5;
					}

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

					PDZDTC = (DDZ1 + DDZ2 + DDZ3) * DT;
					PDZDTC = Math.Max(-0.5, PDZDTC);

					// The change in DZ due to compaction

					DZSNSO[J] *= (1 + PDZDTC);
//					if (double.IsNaN(DZSNSO[0]))
//						throw new Exception("");
				}

				// Pressure of overlying snow

				BURDEN = BURDEN + WX;

			}

			//if(double.IsNaN(DZSNSO[0]))
//				throw new Exception("");
		}
		/// <summary>
		/// SUBROUTINE FRH2O
		/// ----------------------------------------------------------------------
		/// CALCULATE AMOUNT OF SUPERCOOLED LIQUID SOIL WATER CONTENT IF
		/// TEMPERATURE IS BELOW 273.15K (TFRZ).  REQUIRES NEWTON-TYPE ITERATION
		/// TO SOLVE THE NONLINEAR IMPLICIT EQUATION GIVEN IN EQN 17 OF KOREN ET AL
		/// (1999, JGR, VOL 104(D16), 19569-19585).
		/// ----------------------------------------------------------------------
		/// NEW VERSION (JUNE 2001): MUCH FASTER AND MORE ACCURATE NEWTON
		/// ITERATION ACHIEVED BY FIRST TAKING LOG OF EQN CITED ABOVE -- LESS THAN
		/// 4 (TYPICALLY 1 OR 2) ITERATIONS ACHIEVES CONVERGENCE.  ALSO, EXPLICIT
		/// 1-STEP SOLUTION OPTION FOR SPECIAL CASE OF PARAMETER CK=0, WHICH
		/// REDUCES THE ORIGINAL IMPLICIT EQUATION TO A SIMPLER EXPLICIT FORM,
		/// KNOWN AS THE "FLERCHINGER EQN". IMPROVED HANDLING OF SOLUTION IN THE
		/// LIMIT OF FREEZING POINT TEMPERATURE TFRZ.
		/// </summary>
		static double FRH2O(GridCell cell, double TKELV, double SMC, double SH2O)
		{

			// INPUT:

			//   TKELV.........TEMPERATURE (Kelvin)
			//   SMC...........TOTAL SOIL MOISTURE CONTENT (VOLUMETRIC)
			//   SH2O..........LIQUID SOIL MOISTURE CONTENT (VOLUMETRIC)
			//   B.............SOIL TYPE "B" PARAMETER (FROM REDPRM)
			//   PSISAT........SATURATED SOIL MATRIC POTENTIAL (FROM REDPRM)

			// OUTPUT:
			//   FREE..........SUPERCOOLED LIQUID WATER CONTENT [m3/m3]
			// ----------------------------------------------------------------------
			//   REAL, INTENT(IN)     :: SH2O,SMC,TKELV
			//    REAL, INTENT(OUT)    :: FREE
			//    REAL                 :: BX,DENOM,DF,DSWL,FK,SWL,SWLK
			//    INTEGER              :: NLOG,KCOUNT
			//      PARAMETER(CK = 0.0)
			double SWL = 0;
			double CK = 8.0;
			double BLIM = 5.5;
			double ERROR = 0.00005;
			double DICE = 920.0;
			double FREE = 0;
			

			// ----------------------------------------------------------------------
			// LIMITS ON PARAMETER B: B < 5.5  (use parameter BLIM)
			// SIMULATIONS SHOWED IF B > 5.5 UNFROZEN WATER CONTENT IS
			// NON-REALISTICALLY HIGH AT VERY LOW TEMPERATURES.
			// ----------------------------------------------------------------------
			double BX = cell.BEXP;
			// ----------------------------------------------------------------------
			// INITIALIZING ITERATIONS COUNTER AND ITERATIVE SOLUTION FLAG.
			// ----------------------------------------------------------------------

			if (cell.BEXP > BLIM)
				BX = BLIM;
			int NLOG = 0;

			// ----------------------------------------------------------------------
			//  if TEMPERATURE NOT SIGNifICANTLY BELOW FREEZING (TFRZ), SH2O = SMC
			// ----------------------------------------------------------------------
			int KCOUNT = 0;
			if (TKELV > (NoahMP.TFRZ - 1E-3))
				FREE = SMC;
			else {
				// ----------------------------------------------------------------------
				// OPTION 1: ITERATED SOLUTION IN KOREN ET AL, JGR, 1999, EQN 17
				// ----------------------------------------------------------------------
				// INITIAL GUESS FOR SWL (frozen content)
				// ----------------------------------------------------------------------
				if (CK != 0.0) {
					SWL = SMC - SH2O;
					// ----------------------------------------------------------------------
					// KEEP WITHIN BOUNDS.
					// ----------------------------------------------------------------------
					if (SWL > (SMC - 0.02))
						SWL = SMC - 0.02;
					// ----------------------------------------------------------------------
					//  START OF ITERATIONS
					// ----------------------------------------------------------------------
					if (SWL < 0)
						SWL = 0;
					//1001      Continue
					while (true) {
						if (!((NLOG < 10) && (KCOUNT == 0)))
							break;
						
						NLOG = NLOG + 1;
						double	DF = Math.Log((cell.PSISAT * NoahMP.GRAV / NoahMP.HFUS) * Math.Pow(1 + CK * SWL, 2) *
						            Math.Pow(cell.SMCMAX / (SMC - SWL), BX)) - Math.Log(-(TKELV - NoahMP.TFRZ) / TKELV);
						double DENOM = 2 * CK / (1 + CK * SWL) + BX / (SMC - SWL);
						double SWLK = SWL - DF / DENOM;
						// ----------------------------------------------------------------------
						// BOUNDS USEFUL FOR MATHEMATICAL SOLUTION.
						// ----------------------------------------------------------------------
						if (SWLK > (SMC - 0.02))
							SWLK = SMC - 0.02;
						if (SWLK < 0)
							SWLK = 0;
//						if (double.IsNaN(SWLK))
//							throw new Exception();
						// ----------------------------------------------------------------------
						// MATHEMATICAL SOLUTION BOUNDS APPLIED.
						// ----------------------------------------------------------------------
						double DSWL = Math.Abs(SWLK - SWL);
						// if MORE THAN 10 ITERATIONS, USE EXPLICIT METHOD (CK=0 APPROX.)
						// WHEN DSWL LESS OR EQ. ERROR, NO MORE ITERATIONS REQUIRED.
						// ----------------------------------------------------------------------
						SWL = SWLK;
						if (DSWL <= ERROR) {
							KCOUNT += 1;
						}
						
						// ----------------------------------------------------------------------
						//  END OF ITERATIONS
						// ----------------------------------------------------------------------
						// BOUNDS APPLIED WITHIN DO-BLOCK ARE VALID FOR PHYSICAL SOLUTION.
						// ----------------------------------------------------------------------
					}  //end while
					FREE = SMC - SWL;
//					if (double.IsNaN(FREE))
//						throw new Exception();
					
					// ----------------------------------------------------------------------
					// END OPTION 1
					// ----------------------------------------------------------------------
					// ----------------------------------------------------------------------
					// OPTION 2: EXPLICIT SOLUTION FOR FLERCHINGER EQ. i.e. CK=0
					// IN KOREN ET AL., JGR, 1999, EQN 17
					// APPLY PHYSICAL BOUNDS TO FLERCHINGER SOLUTION
					// ----------------------------------------------------------------------
					if (KCOUNT == 0) {
						//          write(message, '("Flerchinger used in NEW version. Iterations=", I6)') NLOG
						//call wrf_message(trim(message))
						double temp = ((NoahMP.HFUS / (NoahMP.GRAV * (-cell.PSISAT))) * ((TKELV - NoahMP.TFRZ) / TKELV));
						double FK = Math.Pow(temp, -1 / BX) * cell.SMCMAX;
						if (FK < 0.02)
							FK = 0.02;
						FREE = Math.Min(FK, SMC);
//						if (double.IsNaN(FREE))
//							throw new Exception();
						
					}
				} //end if
				// ----------------------------------------------------------------------
				// END OPTION 2
				// ----------------------------------------------------------------------
			}
			return FREE;
		}
		public static void CO2FLUX(int NSNOW, int  NSOIL, int VEGTYP, double IGS, double DT,   //in
			FortDoubleArray DZSNSO, FortDoubleArray STC, double PSN, double TROOT, double TV,   //in
			double WROOT, double WSTRES, double FOLN, double LAPM,           //in
			double LAT, double ILOC, double JLOC, double FVEG,            //in
			ref double XLAI, ref double XSAI, ref double LFMASS, ref double RTMASS, ref double STMASS,   //inout
			ref double FASTCP, ref double STBLCP, ref double WOOD,               //inout
			out double     GPP, out double NPP, out double NEE, out double AUTORS, out double HETERS,    //out
			out   double TOTSC, out double TOTLB)                          //out
		{
			// constants
			const double RTOVRC = 2.0E-8;        //original was 2.0e-8   root turnover coefficient [1/s]
			const double RSDRYC = 40.0;          //original was 40.0   degree of drying that reduces soil respiration [-]
			const double RSWOODC = 3.0E-10;       //wood respiration coeficient [1/s]
			const double BF = 0.90;          //original was 0.90   // carbon to roots
			const double WSTRC = 100.0;   //water stress coeficient [-]
			const double LAIMIN = 0.05;  //minimum leaf area index [m2/m2]
			const double XSAMIN = 0.01;  //minimum leaf area index [m2/m2]

			const double SAPM = 3 * 0.001;     // m2/kg -->m2/g    stem area per unit mass (m2/g)
			double LFMSMN = LAIMIN / LAPM;
			const double STMSMN = XSAMIN / SAPM;
			// ---------------------------------------------------------------------------------

			// respiration

			double RF = 1.0;       //respiration reduction factor (<= 1)
			if (Math.Abs(IGS) < 1e-10)
				RF = 0.5;
			//计算维持呼吸
			double FNF = Math.Min(FOLN / Math.Max(1E-06, Driver.vegparams.FOLNMX[VEGTYP - 1]), 1.0); //foliage nitrogen adjustemt to respiration (<= 1)
			double TF = Math.Pow(Driver.vegparams.ARM[VEGTYP - 1], (TV - 298.16) / 10);  //temperature factor
			double RESP = Driver.vegparams.RMF25[VEGTYP - 1] * TF * FNF * XLAI * RF * (1 - WSTRES); // umol/m2/s  leaf respiration [umol/m2/s]
			double RSLEAF = Math.Min(LFMASS / DT, RESP * 12e-6);                         // g/m2/s  leaf maintenance respiration per timestep [g/m2]
			double RSROOT = Driver.vegparams.RMR25[VEGTYP - 1] * (RTMASS * 1E-3) * TF * RF * 12e-6;         // fine root respiration per time step [g/m2]
			double RSSTEM = Driver.vegparams.RMS25[VEGTYP - 1] * (STMASS * 1E-3) * TF * RF * 12e-6;         // g/m2/s	stem respiration [g/m2/s]
			double RSWOOD = RSWOODC * R(TV) * WOOD * Driver.vegparams.WDPOOL[VEGTYP - 1]; //wood respiration [g/m2]

			// carbon assimilation
			// 1 mole -> 12 g carbon or 44 g CO2; 1 umol -> 12.e-6 g carbon;
			//有一个bug,即PSN永远为0
			double CARBFX = PSN * 12e-6;              // umol co2 /m2/ s -> g/m2/s carbon			

			// fraction of carbon into leaf versus nonleaf,即碳进入叶子和非叶子的比例
			//从LAI计算出LEAFPT（叶子中的碳占比）
			double LEAFPT = Math.Exp(0.01 * (1 - Math.Exp(0.75 * XLAI)) * XLAI);  //fraction of carbon allocated to leaves [-]
			if (VEGTYP == Driver.vegparams.EBLFOREST)
				LEAFPT = Math.Exp(0.01 * (1 - Math.Exp(0.50 * XLAI)) * XLAI);

			double NONLEF = 1.0 - LEAFPT;
			double STEMPT = XLAI / 10.0; //茎占LAI中10%的比例
			LEAFPT -= STEMPT;

			//  fraction of carbon into wood versus root
			double WOODF = 0;
			if (WOOD > 0)
				WOODF = (1 - Math.Exp(-BF * (Driver.vegparams.WRRAT[VEGTYP - 1] * RTMASS / WOOD)) / BF) * Driver.vegparams.WDPOOL[VEGTYP - 1];
			else
				WOODF = 0;
			

			double ROOTPT = NONLEF * (1 - WOODF);
			double WOODPT = NONLEF * WOODF;
			
			//衰老和机械死亡: leaf and root turnover per time step
			double LTOVRC = Driver.vegparams.LTOVRC[VEGTYP - 1] / 10; 
			double LFTOVR = LTOVRC * 1E-6 * LFMASS;
			double STTOVR = LTOVRC * 1E-6 * STMASS;
			double RTTOVR = RTOVRC * RTMASS;
			double WDTOVR = 9.5E-10 * WOOD;

			//干冻死亡：seasonal leaf die rate dependent on temp and water stress
			// water stress is set to 1 at permanent wilting point
			double SC = Math.Exp(-0.3 * Math.Max(0, TV - Driver.vegparams.TDLEF[VEGTYP - 1])) * (LFMASS / 120);
			double SD = Math.Exp((WSTRES - 1) * WSTRC);
			//DIELF:death of leaf mass per time step [g/m2]
			double DIELF = LFMASS * 1E-6 * (Driver.vegparams.DILEFW[VEGTYP - 1] * SD + Driver.vegparams.DILEFC[VEGTYP - 1] * SC);
			double DIEST = STMASS * 1E-6 * (Driver.vegparams.DILEFW[VEGTYP - 1] * SD + Driver.vegparams.DILEFC[VEGTYP - 1] * SC);

			//计算生长呼吸 calculate growth respiration for leaf, rtmass and wood:除去用于"维持呼吸"的部分光合作用量，剩下的部分乘以一个因子
			double FRAGR = Driver.vegparams.FRAGR[VEGTYP - 1];
			double GRLEAF = Math.Max(0.0, FRAGR * (LEAFPT * CARBFX - RSLEAF));  //growth respiration rate for leaf [g/m2/s]
			double GRSTEM = Math.Max(0.0, FRAGR * (STEMPT * CARBFX - RSSTEM));  //growth respiration rate for stem [g/m2/s]
			double GRROOT = Math.Max(0.0, FRAGR * (ROOTPT * CARBFX - RSROOT));  //growth respiration rate for root [g/m2/s]
			double GRWOOD = Math.Max(0.0, FRAGR * (WOODPT * CARBFX - RSWOOD));  //growth respiration rate for wood [g/m2/s]

			// Impose lower T limit for photosynthesis
			//减去呼吸消耗
			//ADDNPPLF:leaf assimil after resp. losses removed [g/m2]
			double ADDNPPLF = Math.Max(0, LEAFPT * CARBFX - GRLEAF - RSLEAF);
			double ADDNPPST = Math.Max(0, STEMPT * CARBFX - GRSTEM - RSSTEM);  //stem assimil after resp. losses removed [g/m2]
			//下面是刘永和注释掉的。下面的这几句运行后会使输出的空间分布变得很异样
			if (TV < Driver.vegparams.TMIN[VEGTYP - 1]) {
				ADDNPPLF = 0;
				ADDNPPST = 0;
			}

			// update leaf, root, and wood carbon
			// avoid reducing leaf mass below its minimum value but conserve mass
			
			//干冻死亡消耗
			double LFDEL = (LFMASS - LFMSMN) / DT; //maximum  leaf mass  available to change [g/m2/s]
			double STDEL = (STMASS - STMSMN) / DT;
			//对干冻死亡量要做限制性修正
			DIELF = Math.Min(DIELF, LFDEL + ADDNPPLF - LFTOVR);  
			DIEST = Math.Min(DIEST, STDEL + ADDNPPST - STTOVR);

			// net primary productivities
			double NPPL = Math.Max(ADDNPPLF, -LFDEL);         //leaf net primary productivity [g/m2/s]
			double NPPS = Math.Max(ADDNPPST, -STDEL);         //wood net primary productivity [g/m2/s]
			double NPPR = ROOTPT * CARBFX - RSROOT - GRROOT;  //root net primary productivity [g/m2/s]
			double NPPW = WOODPT * CARBFX - RSWOOD - GRWOOD;  //wood net primary productivity [g/m2/s]

			// masses of plant components
			//减去上面的两项消耗
			LFMASS += (NPPL - LFTOVR - DIELF) * DT;
			STMASS += (NPPS - STTOVR - DIEST) * DT;   // g/m2
			RTMASS += (NPPR - RTTOVR) * DT;
			double delta = (NPPL - LFTOVR - DIELF) * DT;
			//Console.WriteLine("LFMASS=" + LFMASS + " delta=" + delta + " NPPL=" + NPPL);
			if (RTMASS < 0.0) {
				RTTOVR = NPPR;
				RTMASS = 0.0;
			}
			WOOD = (WOOD + (NPPW - WDTOVR) * DT) * Driver.vegparams.WDPOOL[VEGTYP - 1];

//			if(WOOD<0)
//				WOOD=0;          //刘永和加上的
//				throw new Exception("");
			// soil carbon budgets

			FASTCP += (RTTOVR + LFTOVR + STTOVR + WDTOVR + DIELF) * DT;

			//计算土壤碳呼吸
			double FST = Math.Pow(2.0, (STC[1] - 283.16) / 10);
			double FSW = WROOT / (0.20 + WROOT) * 0.23 / (0.23 + WROOT);
			double RSSOIL = FSW * FST * Driver.vegparams.MRP[VEGTYP - 1] * Math.Max(0, FASTCP * 1E-3) * 12E-6;
			double STABLC = 0.1 * RSSOIL;
			FASTCP -= (RSSOIL + STABLC) * DT;
			STBLCP += STABLC * DT;

			//  total carbon flux
			double CFLUX = -CARBFX + RSLEAF + RSROOT + RSWOOD + RSSTEM + RSSOIL + GRLEAF + GRROOT + GRWOOD;// g/m2/s
			

			// for outputs

			GPP = CARBFX;                                             //g/m2/s C
			NPP = NPPL + NPPW + NPPR;                                 //g/m2/s C
			//呼吸消耗总量
			AUTORS = RSROOT + RSWOOD + RSLEAF + GRLEAF + GRROOT + GRWOOD;   //g/m2/s C  //net ecosystem resp. (maintance and growth)
			//g/m2/s C
			HETERS = RSSOIL;                                             //g/m2/s C
			//净生态系统碳交换量
			NEE = (AUTORS + HETERS - GPP) * 44 / 12;                    //g/m2/s CO2
			TOTSC = FASTCP + STBLCP;                                    //g/m2   C
			TOTLB = LFMASS + RTMASS + WOOD;                             //g/m2   C
			//Console.WriteLine("NPP="+NPP+"	GPP="+GPP+"	WOOD="+WOOD);

			// leaf area index and stem area index
			XLAI = Math.Max(LFMASS * LAPM, LAIMIN);
			XSAI = Math.Max(STMASS * SAPM, XSAMIN);
//			if(PSN>0)
//				Console.WriteLine("PSN="+PSN+" XLAI="+XLAI);
			
		}
		/// <summary>
		/// CO2通量(参照WRF-Hydro5.2)，该过程已验证过一次，结果一致，只是有一个指数次方造成的较小误差
		/// </summary>
		/// <param name="NSNOW"></param>
		/// <param name="NSOIL"></param>
		/// <param name="VEGTYP"></param>
		/// <param name="IGS">growing season index (0=off, 1=on)</param>
		/// <param name="DT"></param>
		/// <param name="DZSNSO"></param>
		/// <param name="STC"></param>
		/// <param name="PSN">total photosyn. (umolco2/m2/s) [+]，总光合作用</param>
		/// <param name="TROOT"></param>
		/// <param name="TV">leaf temperature (k)</param>
		/// <param name="WROOT">root zone soil water</param>
		/// <param name="WSTRES">soil water stress</param>
		/// <param name="FOLN">foliage nitrogen (%)</param>
		/// <param name="LAPM">leaf area per unit mass [m2/g]，单位质量的叶片面积，与具体植被有关</param>
		/// <param name="LAT"></param>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="FVEG"></param>
		/// <param name="XLAI">叶面积指数</param>
		/// <param name="XSAI">茎面积指数</param>
		/// <param name="LFMASS"></param>
		/// <param name="RTMASS"></param>
		/// <param name="STMASS"></param>
		/// <param name="FASTCP">short lived carbon [g/m2]</param>
		/// <param name="STBLCP">stable carbon pool [g/m2]</param>
		/// <param name="WOOD"></param>
		/// <param name="GPP">总初级生产力</param>
		/// <param name="NPP">net primary productivity [g/m2] 净初级生产力</param>
		/// <param name="NEE">net ecosystem exchange (autors+heters-gpp) [g/m2/s CO2] 净生态系统碳交换量</param>
		/// <param name="AUTORS">net ecosystem resp. (maintance and growth)</param>
		/// <param name="HETERS">organic respiration</param>
		/// <param name="TOTSC">total soil carbon (g/m2)</param>
		/// <param name="TOTLB">total living carbon (g/m2)</param>
		public static void CO2FLUX_52(int NSNOW, int  NSOIL, int VEGTYP, double IGS, double DT,   //in
			FortDoubleArray DZSNSO, FortDoubleArray STC, double PSN, double TROOT, double TV,   //in
			double WROOT, double WSTRES, double FOLN, double LAPM,           //in
			double LAT, double ILOC, double JLOC, double FVEG,            //in
			ref double XLAI, ref double XSAI, ref double LFMASS, ref double RTMASS, ref double STMASS,   //inout
			ref double FASTCP, ref double STBLCP, ref double WOOD,               //inout
			out double GPP, out double NPP, out double NEE, out double AUTORS, out double HETERS,    //out
			out   double TOTSC, out double TOTLB)                          //out
		{
			// constants
			const double RTOVRC = 2.0E-8;        //original was 2.0e-8   root turnover coefficient [1/s]
			const double RSDRYC = 40.0;          //original was 40.0   degree of drying that reduces soil respiration [-]
			const double RSWOODC = 3.0E-10;       //wood respiration coeficient [1/s]
			const double BF = 0.90;          //original was 0.90   // carbon to roots
			const double WSTRC = 100.0;   //water stress coeficient [-]
			const double LAIMIN = 0.05;  //minimum leaf area index [m2/m2]
			const double XSAMIN = 0.05;  //minimum leaf area index [m2/m2]
			const double SAPM = 3 * 0.001;     // m2/kg -->m2/g    stem area per unit mass (m2/g)
			double LFMSMN = LAIMIN / LAPM;     //最小叶子量。LAPM为单位物质的分布面积
			const double STMSMN = XSAMIN / SAPM;
			// ---------------------------------------------------------------------------------

			double RF = 1.0;       //respiration reduction factor (<= 1)
			if (Math.Abs(IGS) < 1e-10)
				RF = 0.5;
			//叶子氮系数：foliage nitrogen adjustemt to respiration (<= 1), 0.67
			double FNF = Math.Min(FOLN / Math.Max(1E-06, Driver.vegparams.FOLNMX[VEGTYP - 1]), 1.0); 
			double TF = Math.Pow(Driver.vegparams.ARM[VEGTYP - 1], (TV - 298.16) / 10);  //temperature factor, 即alpha
			//下面就是梁晶晶论文中的公式。WSTRES为土壤水胁迫比，等于0时最优，1时最差
			double RESP = Driver.vegparams.RMF25[VEGTYP - 1] * TF * FNF * XLAI * RF * (1 - WSTRES); // umol/m2/s  leaf respiration [umol/m2/s]。
			//参WRF-Hydro5.2进行了修改
			double RSLEAF = Math.Min((LFMASS - LFMSMN) / DT, RESP * 12e-6);                         // g/m2/s  leaf maintenance respiration per timestep [g/m2]
			double RSROOT = Driver.vegparams.RMR25[VEGTYP - 1] * (RTMASS * 1E-3) * TF * RF * 12e-6;         // fine root respiration per time step [g/m2]
			double RSSTEM = Driver.vegparams.RMS25[VEGTYP - 1] * (STMASS * 1E-3) * TF * RF * 12e-6;         // g/m2/s	stem respiration [g/m2/s]
			double RSWOOD = RSWOODC * R(TV) * WOOD * Driver.vegparams.WDPOOL[VEGTYP - 1]; //wood respiration [g/m2]

			// carbon assimilation
			// 1 mole -> 12 g carbon or 44 g CO2; 1 umol -> 12.e-6 g carbon;
			
			double CARBFX = PSN * 12e-6;              // umol co2 /m2/ s -> g/m2/s carbon			

			// fraction of carbon into leaf versus nonleaf,即碳进入叶子和非叶子的比例
			//从LAI计算出LEAFPT（叶子中的碳占比）
			double LEAFPT = Math.Exp(0.01 * (1 - Math.Exp(0.75 * XLAI)) * XLAI);  //fraction of carbon allocated to leaves [-]
			if (VEGTYP == Driver.vegparams.EBLFOREST)
				LEAFPT = Math.Exp(0.01 * (1 - Math.Exp(0.50 * XLAI)) * XLAI);

			double NONLEF = 1.0 - LEAFPT;
			double STEMPT = XLAI * 0.1 * LEAFPT; //茎占LAI中10%的比例，参WRF-Hydro5.2进行了修改
			LEAFPT -= STEMPT;

			//木质比： fraction of carbon into wood versus root
			double WOODF = 0;
			if (WOOD > 0)
				WOODF = (1 - Math.Exp(-BF * (Driver.vegparams.WRRAT[VEGTYP - 1] * RTMASS / WOOD)) / BF) * Driver.vegparams.WDPOOL[VEGTYP - 1];
			else
				WOODF = Driver.vegparams.WDPOOL[VEGTYP - 1];
			

			//非叶子又分为根系和木质
			//根系与木质的占比
			double ROOTPT = NONLEF * (1 - WOODF);
			double WOODPT = NONLEF * WOODF;
			
			//叶子与根系的倒伏 leaf and root turnover per time step]
			double LTOVRC = Driver.vegparams.LTOVRC[VEGTYP - 1] / 2; //这个参数偏大可能会造成异常，是一个非常敏感的参数
			double LFTOVR = LTOVRC * 5E-7 * LFMASS;  //参照WRF-Hydro5.2进行了修改。
			double STTOVR = LTOVRC * 5E-7 * STMASS;
			double RTTOVR = RTOVRC * RTMASS;
			double WDTOVR = 9.5E-10 * WOOD;

			// seasonal leaf die rate dependent on temp and water stress
			// water stress is set to 1 at permanent wilting point
			double SC = Math.Exp(-0.3 * Math.Max(0, TV - Driver.vegparams.TDLEF[VEGTYP - 1])) * (LFMASS / 120);
			double SD = Math.Exp((WSTRES - 1) * WSTRC);
			//DIELF:death of leaf mass per time step [g/m2]
			double DIELF = LFMASS * 1E-6 * (Driver.vegparams.DILEFW[VEGTYP - 1] * SD + Driver.vegparams.DILEFC[VEGTYP - 1] * SC);
			double DIEST = STMASS * 1E-6 * (Driver.vegparams.DILEFW[VEGTYP - 1] * SD + Driver.vegparams.DILEFC[VEGTYP - 1] * SC);

			//计算生长呼吸 calculate growth respiration for leaf, rtmass and wood:除去用于"维持呼吸"的部分光合作用量，剩下的部分乘以一个因子
			double FRAGR = Driver.vegparams.FRAGR[VEGTYP - 1];
			double GRLEAF = Math.Max(0.0, FRAGR * (LEAFPT * CARBFX - RSLEAF));  //growth respiration rate for leaf [g/m2/s]
			double GRSTEM = Math.Max(0.0, FRAGR * (STEMPT * CARBFX - RSSTEM));  //growth respiration rate for stem [g/m2/s]
			double GRROOT = Math.Max(0.0, FRAGR * (ROOTPT * CARBFX - RSROOT));  //growth respiration rate for root [g/m2/s]
			double GRWOOD = Math.Max(0.0, FRAGR * (WOODPT * CARBFX - RSWOOD));  //growth respiration rate for wood [g/m2/s]

			// Impose lower T limit for photosynthesis
			//减去呼吸消耗
			//ADDNPPLF:leaf assimil after resp. losses removed [g/m2]
			double ADDNPPLF = Math.Max(0, LEAFPT * CARBFX - GRLEAF - RSLEAF);
			double ADDNPPST = Math.Max(0, STEMPT * CARBFX - GRSTEM - RSSTEM);  //stem assimil after resp. losses removed [g/m2]
			
			if (TV < Driver.vegparams.TMIN[VEGTYP - 1]) {
				ADDNPPLF = 0;
				ADDNPPST = 0;
			}

			// update leaf, root, and wood carbon
			// avoid reducing leaf mass below its minimum value but conserve mass
			
			//温度和水供应导致的植被死亡,不能超过最大死亡量
			double LFDEL = (LFMASS - LFMSMN) / DT; //maximum  leaf mass  available to change [g/m2/s]
			double STDEL = (STMASS - STMSMN) / DT;
			//对干冻死亡量要做限制性修正
			DIELF = Math.Min(DIELF, LFDEL + ADDNPPLF - LFTOVR);
			DIEST = Math.Min(DIEST, STDEL + ADDNPPST - STTOVR);

			// net primary productivities. NPP=光合作用-(maintenance呼吸+生长呼吸),因为植物在夜间也会释放出CO2
			double NPPL = Math.Max(ADDNPPLF, -LFDEL);         //leaf net primary productivity [g/m2/s]
			double NPPS = Math.Max(ADDNPPST, -STDEL);         //stem net primary productivity [g/m2/s]
			double NPPR = ROOTPT * CARBFX - RSROOT - GRROOT;  //root net primary productivity [g/m2/s]
			double NPPW = WOODPT * CARBFX - RSWOOD - GRWOOD;  //wood net primary productivity [g/m2/s]

			// masses of plant components
			//再减去两项消耗：冷干死亡(DIELF)+衰老死亡(LFTOVR)
			LFMASS += (NPPL - LFTOVR - DIELF) * DT;
			STMASS += (NPPS - STTOVR - DIEST) * DT;   // g/m2
			RTMASS += (NPPR - RTTOVR) * DT;
			double delta = (NPPL - LFTOVR - DIELF) * DT;
			//Console.WriteLine("LFMASS=" + LFMASS + " delta=" + delta + " NPPL=" + NPPL);
			if (RTMASS < 0.0) {
				RTTOVR = NPPR;
				RTMASS = 0.0;
			}
			WOOD = (WOOD + (NPPW - WDTOVR) * DT) * Driver.vegparams.WDPOOL[VEGTYP - 1]; //WDPOOL是0或1.


			//FASTCP加上所有的死亡量， soil carbon budgets
			FASTCP += (RTTOVR + LFTOVR + STTOVR + WDTOVR + DIELF) * DT;

			double FST = Math.Pow(2.0, (STC[1] - 283.16) / 10);
			double FSW = WROOT / (0.20 + WROOT) * 0.23 / (0.23 + WROOT);
			double RSSOIL = FSW * FST * Driver.vegparams.MRP[VEGTYP - 1] * Math.Max(0, FASTCP * 1E-3) * 12E-6;
			double STABLC = 0.1 * RSSOIL;
			FASTCP -= (RSSOIL + STABLC) * DT;
			STBLCP += STABLC * DT;

			//  total carbon flux, 参照WRF-Hydro5.2进行了修改。
			double CFLUX = -CARBFX + RSLEAF + RSROOT + RSWOOD + RSSTEM
			               + RSSOIL * 0.9 + GRLEAF + GRROOT + GRWOOD + GRSTEM;// g/m2/s  
			
			// for outputs
			GPP = CARBFX;                                             //g/m2/s C
			NPP = NPPL + NPPW + NPPR + NPPS;  //g/m2/s C
			//总呼吸消耗量（维持呼吸+生长呼吸），net ecosystem resp. (maintance and growth)
			AUTORS = RSROOT + RSWOOD + RSLEAF + GRSTEM + GRLEAF + GRROOT + GRWOOD + GRSTEM;   //参照WRF-Hydro5.2进行了修改,g/m2/s C  
			
			//organic respiration
			HETERS = RSSOIL * 0.9;    //参照WRF-Hydro5.2进行了修改, g/m2/s C
			//net ecosystem exchange (autors+heters-gpp)
			NEE = (AUTORS + HETERS - GPP) * 44 / 12;                    //g/m2/s CO2
			//total soil carbon (g/m2)
			TOTSC = FASTCP + STBLCP;                                    //g/m2   C
			//total living carbon (g/m2)
			TOTLB = LFMASS + RTMASS + STMASS + WOOD;  //参照WRF-Hydro5.2进行了修改, g/m2   C			

			// leaf area index and stem area index
			XLAI = Math.Max(LFMASS * LAPM, LAIMIN);
			XSAI = Math.Max(STMASS * SAPM, XSAMIN);
//			Console.WriteLine("PSN="+PSN+" XLAI="+XLAI+" XSAI="+XSAI);

		}
		public static void CSNOW(int ISNOW, int NSNOW, int NSOIL, FortDoubleArray SNICE, FortDoubleArray SNLIQ, FortDoubleArray DZSNSO, FortDoubleArray TKSNO, FortDoubleArray CVSNO, FortDoubleArray SNICEV, FortDoubleArray SNLIQV, FortDoubleArray EPORE)
		{
			// --------------------------------------------------------------------------------------------------
			// Snow bulk density,volumetric capacity, and thermal conductivity
			//---------------------------------------------------------------------------------------------------
			// inputs

			//  INTEGER,                          INTENT(IN) :: ISNOW  //number of snow layers (-)
			//  INTEGER                        ,  INTENT(IN) :: NSNOW  //maximum no. of snow layers
			//  INTEGER                        ,  INTENT(IN) :: NSOIL  //number of soil layers
			//  REAL, DIMENSION(-NSNOW+1:    0),  INTENT(IN) :: SNICE  //snow ice mass (kg/m2)
			//  REAL, DIMENSION(-NSNOW+1:    0),  INTENT(IN) :: SNLIQ  //snow liq mass (kg/m2)
			//  REAL, DIMENSION(-NSNOW+1:NSOIL),  INTENT(IN) :: DZSNSO //snow/soil layer thickness [m]

			// outputs

			//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: CVSNO  //volumetric specific heat (j/m3/k)
			//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: TKSNO  //thermal conductivity (w/m/k)
			//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: SNICEV //partial volume of ice [m3/m3]
			//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: SNLIQV //partial volume of liquid water [m3/m3]
			//  REAL, DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: EPORE  //effective porosity [m3/m3]

			// locals

			//  INTEGER :: IZ
			//  REAL, DIMENSION(-NSNOW+1:    0) :: BDSNOI  //bulk density of snow(kg/m3)
			
			// thermal capacity of snow
			FortDoubleArray BDSNOI = new FortDoubleArray(1 - NSNOW, 0);
			
			for (int IZ = ISNOW + 1; IZ <= 0; IZ++) {
				double temp = Math.Max(DZSNSO[IZ], 1e-6);
				SNICEV[IZ] = Math.Min(1, SNICE[IZ] / (temp * NoahMP.DENICE));
				EPORE[IZ] = 1 - SNICEV[IZ];
				SNLIQV[IZ] = Math.Min(EPORE[IZ], SNLIQ[IZ] / (temp * NoahMP.DENH2O));
			}

			for (int IZ = ISNOW + 1; IZ <= 0; IZ++) {
				double temp = Math.Max(DZSNSO[IZ], 1e-6);
				BDSNOI[IZ] = (SNICE[IZ] + SNLIQ[IZ]) / temp;
				CVSNO[IZ] = NoahMP.CICE * SNICEV[IZ] + NoahMP.CWAT * SNLIQV[IZ];
				//      CVSNO[IZ] = 0.525E06;                          // constant
			}

			// thermal conductivity of snow

			for (int IZ = ISNOW + 1; IZ <= 0; IZ++) {
				TKSNO[IZ] = 3.2217E-6 * Math.Pow(BDSNOI[IZ], 2);           // Stieglitz(yen,1965)
				//    TKSNO[IZ] = 2E-2+2.5E-6*BDSNOI[IZ]*BDSNOI[IZ]   // Anderson, 1976
				//    TKSNO[IZ] = 0.35                                // constant
				//    TKSNO[IZ] = 2.576E-6*BDSNOI[IZ]**2. + 0.074    // Verseghy (1991)
				//    TKSNO[IZ] = 2.22*(BDSNOI[IZ]/1000)**1.88      // Douvill(Yen, 1981)
				//if(double.IsNaN( TKSNO[IZ]))
//	throw new Exception("");
			}
		}
		public static void TWOSTREAM(int IB, int IC, int VEGTYP, double COSZ, double VAI,
			double FWET, double T, double[] ALBGRD, double[] ALBGRI, double[] RHO,
			double[] TAU, double FVEG, int IST, int ILOC, int JLOC,
			double[] FAB, double[] FRE, double[] FTD, double[] FTI, out double GDIR, // //)   //out
			double[] FREV, double[] FREG, out double BGAP, out double WGAP)
		{
			// input

			//   INTEGER,              INTENT(IN)  :: ILOC    //grid index
			//   INTEGER,              INTENT(IN)  :: JLOC    //grid index
			//   INTEGER,              INTENT(IN)  :: IST     //surface type
			//   INTEGER,              INTENT(IN)  :: IB      //waveband number
			//   INTEGER,              INTENT(IN)  :: IC      //0=unit incoming direct; 1=unit incoming diffuse
			//   INTEGER,              INTENT(IN)  :: VEGTYP  //vegetation type
//
			//   REAL,                 INTENT(IN)  :: COSZ    //cosine of direct zenith angle (0-1)
			//   REAL,                 INTENT(IN)  :: VAI     //one-sided leaf+stem area index (m2/m2)
			//   REAL,                 INTENT(IN)  :: FWET    //fraction of lai, sai that is wetted (-)
			//   REAL,                 INTENT(IN)  :: T       //surface temperature (k)
//
			//   REAL, DIMENSION(1:2), INTENT(IN)  :: ALBGRD  //direct  albedo of underlying surface (-)
			//   REAL, DIMENSION(1:2), INTENT(IN)  :: ALBGRI  //diffuse albedo of underlying surface (-)
			//   REAL, DIMENSION(1:2), INTENT(IN)  :: RHO     //leaf+stem reflectance
			//   REAL, DIMENSION(1:2), INTENT(IN)  :: TAU     //leaf+stem transmittance
			//   REAL,                 INTENT(IN)  :: FVEG    //green vegetation fraction [0.0-1.0]
//
			// output
//
			//   REAL, DIMENSION(1:2), INTENT(OUT) :: FAB     //flux abs by veg layer (per unit incoming flux)
			//   REAL, DIMENSION(1:2), INTENT(OUT) :: FRE     //flux refl above veg layer (per unit incoming flux)
			//   REAL, DIMENSION(1:2), INTENT(OUT) :: FTD     //down dir flux below veg layer (per unit in flux)
			//   REAL, DIMENSION(1:2), INTENT(OUT) :: FTI     //down dif flux below veg layer (per unit in flux)
			//   REAL,                 INTENT(OUT) :: GDIR    //projected leaf+stem area in solar direction
			//   REAL, DIMENSION(1:2), INTENT(OUT) :: FREV    //flux reflected by veg layer   (per unit incoming flux)
			//   REAL, DIMENSION(1:2), INTENT(OUT) :: FREG    //flux reflected by ground (per unit incoming flux)


			// local
			double OMEGA = 0;     //fraction of intercepted radiation that is scattered
			//   REAL                              :: OMEGAL  //omega for leaves
			//   REAL                              :: BETAI   //upscatter parameter for diffuse radiation
			//   REAL                              :: BETAIL  //betai for leaves
			double BETAD = 0;   //upscatter parameter for direct beam radiation
			//   REAL                              :: BETADL  //betad for leaves
			//   REAL                              :: EXT     //optical depth of direct beam per unit leaf area
			//   REAL                              :: AVMU    //average diffuse optical depth
//
			//   REAL                              :: COSZI   //0.001 <= cosz <= 1.000
			//   REAL                              :: ASU     //single scattering albedo
			//   REAL                              :: CHIL    // -0.4 <= xl <= 0.6
//
			//double TMP0,TMP1,TMP2,TMP3,TMP4,TMP5,TMP6,TMP7,TMP8,TMP9;
			//   REAL                              :: P1,P2,P3,P4,S1,S2,U1,U2,U3
			//   REAL                              :: B,C,D,D1,D2,F,H,H1,H2,H3,H4,H5,H6,H7,H8,H9,H10
			//   REAL                              :: PHI1,PHI2,SIGMA
			double FTDS, FTIS, FRES;
			double DENFVEG = 0;
			//   REAL                              :: VAI_SPREAD
			//jref:start
			double FREVEG, FREBAR, FTDVEG, FTIVEG, FTDBAR, FTIBAR;
			//   REAL                              :: THETAZ
			//  variables for the modified two-stream scheme
			//  Niu and Yang (2004), JGR
//
			double PAI = Math.PI;
			double HD = 0; //     //crown depth (m)
			//   REAL :: BB       //vertical crown radius (m)
			//   REAL :: THETAP   //angle conversion from SZA
			//   REAL :: FA       //foliage volume density (m-1)
			//   REAL :: NEWVAI   //effective LSAI (-)
//
			//   REAL,INTENT(INOUT) :: BGAP     //between canopy gap fraction for beam (-)
			//   REAL,INTENT(INOUT) :: WGAP     //within canopy gap fraction for beam (-)
//
			double KOPEN = 0;   // //gap fraction for diffue light (-)
			double GAP = 0;  //      //total gap fraction for beam ( <=1-shafac )
			// compute within and between gaps
			WGAP = 0;
			BGAP = 0;
			double VAI_SPREAD = VAI;
			if (VAI == 0.0) {
				GAP = 1.0;
				KOPEN = 1.0;
			} else {
				if (NoahMP.OPT_RAD == 1) { //此选项下需要再验证一下
					DENFVEG = -Math.Log(Math.Max(1.0 - FVEG, 0.01)) / (Math.PI * Math.Pow(Driver.vegparams.RC[VEGTYP - 1], 2));
					HD = Driver.vegparams.HVT[VEGTYP - 1] - Driver.vegparams.HVB[VEGTYP - 1];
					HD = Math.Max(0.01, HD);
					double	BB = 0.5 * HD;
					double	THETAP = Math.Atan(BB / Driver.vegparams.RC[VEGTYP - 1] * Math.Tan(Math.Acos(Math.Max(0.01, COSZ))));
					// BGAP    = Math.Exp(-DEN[VEGTYP-1] * PAI * RC[VEGTYP-1]**2/COS(THETAP) );
					BGAP = Math.Exp(-DENFVEG * PAI * Math.Pow(Driver.vegparams.RC[VEGTYP - 1], 2) / Math.Cos(THETAP));
					double FA = VAI / (1.33 * PAI * Math.Pow(Driver.vegparams.RC[VEGTYP - 1], 3.0) * (BB / Driver.vegparams.RC[VEGTYP - 1]) * DENFVEG);
					double NEWVAI = HD * FA;
					WGAP = (1.0 - BGAP) * Math.Exp(-0.5 * NEWVAI / COSZ);
//					if (double.IsNaN(WGAP))
//						throw new Exception();
					GAP = Math.Min(1.0 - FVEG, BGAP + WGAP);

					KOPEN = 0.05;
				}

				if (NoahMP.OPT_RAD == 2) {
					GAP = 0.0;
					KOPEN = 0.0;
				}

				if (NoahMP.OPT_RAD == 3) {
					GAP = 1.0 - FVEG;
					KOPEN = 1.0 - FVEG;
				}
			}
			// calculate two-stream parameters OMEGA, BETAD, BETAI, AVMU, GDIR, EXT.
			// OMEGA, BETAD, BETAI are adjusted for snow. values for OMEGA*BETAD
			// and OMEGA*BETAI are calculated and  divided by the new OMEGA
			// because the product OMEGA*BETAI, OMEGA*BETAD is used in solution.
			// also, the transmittances and reflectances (TAU, RHO) are linear
			// weights of leaf and stem values.

			double COSZI = Math.Max(0.001, COSZ);
			double CHIL = Math.Min(Math.Max(Driver.vegparams.XL[VEGTYP - 1], -0.4), 0.6);
			if (Math.Abs(CHIL) <= 0.01)
				CHIL = 0.01;
			double PHI1 = 0.5 - 0.633 * CHIL - 0.330 * CHIL * CHIL;
			double PHI2 = 0.877 * (1 - 2 * PHI1);
			GDIR = PHI1 + PHI2 * COSZI;
			double EXT = GDIR / COSZI;
			double AVMU = (1 - PHI1 / PHI2 * Math.Log((PHI1 + PHI2) / PHI1)) / PHI2;
			double OMEGAL = RHO[IB] + TAU[IB];
			double TMP0 = GDIR + PHI2 * COSZI;
			double TMP1 = PHI1 * COSZI;
			double TMP2 = 0;
			double ASU = 0.5 * OMEGAL * GDIR / TMP0 * (1 - TMP1 / TMP0 * Math.Log((TMP1 + TMP0) / TMP1));
			double BETADL = (1 + AVMU * EXT) / (OMEGAL * AVMU * EXT) * ASU;
			double BETAIL = 0.5 * (RHO[IB] + TAU[IB] + (RHO[IB] - TAU[IB]) * Math.Pow((1 + CHIL) / 2, 2)) / OMEGAL;
			// adjust omega, betad, and betai for intercepted snow

			if (T > NoahMP.TFRZ) {                                 //no snow
				TMP0 = OMEGAL;
				TMP1 = BETADL;
				TMP2 = BETAIL;
			} else {
				TMP0 = (1 - FWET) * OMEGAL + FWET * RAD_PARAMS.OMEGAS[IB];
				TMP1 = ((1 - FWET) * OMEGAL * BETADL + FWET * RAD_PARAMS.OMEGAS[IB] * RAD_PARAMS.BETADS) / TMP0;
				TMP2 = ((1 - FWET) * OMEGAL * BETAIL + FWET * RAD_PARAMS.OMEGAS[IB] * RAD_PARAMS.BETAIS) / TMP0;
			}

			OMEGA = TMP0;
			BETAD = TMP1;
			double BETAI = TMP2;
			// absorbed, reflected, transmitted fluxes per unit incoming radiation

			double B = 1 - OMEGA + OMEGA * BETAI;
			double C = OMEGA * BETAI;
			TMP0 = AVMU * EXT;
			double D = TMP0 * OMEGA * BETAD;
			double F = TMP0 * OMEGA * (1 - BETAD);
			TMP1 = B * B - C * C;
			double H = Math.Sqrt(TMP1) / AVMU;
			double SIGMA = TMP0 * TMP0 - TMP1;
			//if (Math.Abs(SIGMA) < 1e-6 ) SIGMA = SIGN(1e-6,SIGMA);
			double P1 = B + AVMU * H;
			double P2 = B - AVMU * H;
			double P3 = B + TMP0;
			double P4 = B - TMP0;
			double S1 = Math.Exp(-H * VAI);
			double S2 = Math.Exp(-EXT * VAI);
			double U1 = 0;
			double U2 = 0;
			double U3 = 0;
			if (IC == 0) {
				U1 = B - C / ALBGRD[IB];
				U2 = B - C * ALBGRD[IB];
				U3 = F + C * ALBGRD[IB];
			} else {
				U1 = B - C / ALBGRI[IB];
				U2 = B - C * ALBGRI[IB];
				U3 = F + C * ALBGRI[IB];
			}
			TMP2 = U1 - AVMU * H;
			double TMP3 = U1 + AVMU * H;
			double D1 = P1 * TMP2 / S1 - P2 * TMP3 * S1;
			double TMP4 = U2 + AVMU * H;
			double TMP5 = U2 - AVMU * H;
			double D2 = TMP4 / S1 - TMP5 * S1;
			double H1 = -D * P4 - C * F;
			double TMP6 = D - H1 * P3 / SIGMA;
			double TMP7 = (D - C - H1 / SIGMA * (U1 + TMP0)) * S2;
			double H2 = (TMP6 * TMP2 / S1 - P2 * TMP7) / D1;
			double H3 = -(TMP6 * TMP3 * S1 - P1 * TMP7) / D1;
			double H4 = -F * P3 - C * D;
			double TMP8 = H4 / SIGMA;
			double TMP9 = (U3 - TMP8 * (U2 - TMP0)) * S2;
			double H5 = -(TMP8 * TMP4 / S1 + TMP9) / D2;
			double H6 = (TMP8 * TMP5 * S1 + TMP9) / D2;
			double H7 = (C * TMP2) / (D1 * S1);
			double H8 = (-C * TMP3 * S1) / D1;
			double H9 = TMP4 / (D2 * S1);
			double H10 = (-TMP5 * S1) / D2;

			// downward direct and diffuse fluxes below vegetation
			// Niu and Yang (2004), JGR.

			if (IC == 0) {
				FTDS = S2 * (1.0 - GAP) + GAP;
				FTIS = (H4 * S2 / SIGMA + H5 * S1 + H6 / S1) * (1.0 - GAP);
			} else {
				FTDS = 0;
				FTIS = (H9 * S1 + H10 / S1) * (1.0 - KOPEN) + KOPEN;
			}
			if (double.IsNaN(FTDS))
				throw new Exception();
			FTD[IB] = FTDS;
			FTI[IB] = FTIS;
			// flux reflected by the surface (veg. and ground)

			if (IC == 0) {
				FRES = (H1 / SIGMA + H2 + H3) * (1.0 - GAP) + ALBGRD[IB] * GAP;
				FREVEG = (H1 / SIGMA + H2 + H3) * (1.0 - GAP);
				FREBAR = ALBGRD[IB] * GAP;                   //jref - separate veg. and ground reflection
			} else {
				FRES = (H7 + H8) * (1.0 - KOPEN) + ALBGRI[IB] * KOPEN;
				FREVEG = (H7 + H8) * (1.0 - KOPEN) + ALBGRI[IB] * KOPEN;
				FREBAR = 0;                                //jref - separate veg. and ground reflection
			}
			FRE[IB] = FRES;

			FREV[IB] = FREVEG;
			FREG[IB] = FREBAR;

			// flux absorbed by vegetation

			FAB[IB] = 1 - FRE[IB] - (1 - ALBGRD[IB]) * FTD[IB] - (1 - ALBGRI[IB]) * FTI[IB];
//			if (double.IsNaN(FAB[IB]))
//				throw new Exception("");
			//if(iloc == 1&&jloc ==  2)
			//  write(*,'(a7,2i2,5(a6,f8.4),2(a9,f8.4))') "ib,ic: ",ib,ic," GAP: ",GAP," FTD: ",FTD[IB]," FTI: ",FTI[IB]," FRE: ",
			//         FRE[IB]," FAB: ",FAB[IB]," ALBGRD: ",ALBGRD[IB]," ALBGRI: ",ALBGRI[IB]
			//end if


		}
		public static double TDFCND(GridCell cell, double SMC, double SH2O)
		{
			// --------------------------------------------------------------------------------------------------
			// Calculate thermal diffusivity and conductivity of the soil.
			// Peters-Lidard approach (Peters-Lidard et al., 1998)
			// --------------------------------------------------------------------------------------------------
			// Code history:
			// June 2001 changes: frozen soil condition.
			// --------------------------------------------------------------------------------------------------
			//    IMPLICIT NONE
			//    REAL, INTENT(IN)       :: SMC    // total soil water
			//    REAL, INTENT(IN)       :: SH2O   // liq. soil water
			//    REAL, INTENT(OUT)      :: DF     // thermal diffusivity

			// local variables
			//    REAL  :: AKE
			//    REAL  :: GAMMD
			//    REAL  :: THKDRY
			//    REAL  :: THKO     // thermal conductivity for other soil components
			//    REAL  :: THKQTZ   // thermal conductivity for quartz
			//    REAL  :: THKSAT   // 
			//    REAL  :: THKS     // thermal conductivity for the solids
			//    REAL  :: THKW     // water thermal conductivity
			//    REAL  :: SATRATIO
			//    REAL  :: XU
			//    REAL  :: XUNFROZ
			// --------------------------------------------------------------------------------------------------
			// We now get quartz as an input argument (set in routine redprm):
			//      DATA QUARTZ /0.82, 0.10, 0.25, 0.60, 0.52,
			//                  0.35, 0.60, 0.40, 0.82/
			// --------------------------------------------------------------------------------------------------
			// If the soil has any moisture content compute a partial sum/product
			// otherwise use a constant value which works well with most soils
			// --------------------------------------------------------------------------------------------------
			//  QUARTZ ....QUARTZ CONTENT (SOIL TYPE DEPENDENT)
			// --------------------------------------------------------------------------------------------------
			// USE AS IN PETERS-LIDARD, 1998 (MODIF. FROM JOHANSEN, 1975).

			//                                  PABLO GRUNMANN, 08/17/98
			// Refs.:
			//      Farouki, O.T.,1986: Thermal properties of soils. Series on Rock
			//              and Soil Mechanics, Vol. 11, Trans Tech, 136 pp.
			//      Johansen, O., 1975: Thermal conductivity of soils. PH.D. Thesis,
			//              University of Trondheim,
			//      Peters-Lidard, C. D., et al., 1998: The effect of soil thermal
			//              conductivity parameterization on surface energy fluxes
			//              and temperatures. Journal of The Atmospheric Sciences,
			//              Vol. 55, pp. 1209-1224.
			// --------------------------------------------------------------------------------------------------
			// NEEDS PARAMETERS
			// POROSITY(SOIL TYPE):
			//      POROS = SMCMAX
			// SATURATION RATIO:
			// PARAMETERS  W/(M.K)
			double SATRATIO = SMC / cell.SMCMAX;
			double THKW = 0.57;
			//      if (QUARTZ <= 0.2) THKO = 3.0
			double THKO = 2.0;
			// SOLIDS' CONDUCTIVITY
			// QUARTZ' CONDUCTIVITY
			double THKQTZ = 7.7;

			// UNFROZEN FRACTION (FROM 1.0, i.e., 100%LIQUID, TO 0. (100% FROZEN))
			double THKS = Math.Pow(THKQTZ, cell.QUARTZ) * Math.Pow(THKO, (1 - cell.QUARTZ));

//			if(SMC==0)
//				throw new Exception("");
			// UNFROZEN VOLUME FOR SATURATION (POROSITY*XUNFROZ)
			double XUNFROZ = 0.001;
			if (SMC != 0)
				XUNFROZ =	SH2O / SMC;
			// SATURATED THERMAL CONDUCTIVITY
			double XU = XUNFROZ * cell.SMCMAX;

			// DRY DENSITY IN KG/M3
			double THKSAT = Math.Pow(THKS, (1 - cell.SMCMAX)) * Math.Pow(NoahMP.TKICE, (cell.SMCMAX - XU)) * Math.Pow(THKW, XU);

			// DRY THERMAL CONDUCTIVITY IN W.M-1.K-1
			double GAMMD = (1 - cell.SMCMAX) * 2700;

			double THKDRY = (0.135 * GAMMD + 64.7) / (2700 - 0.947 * GAMMD);
			// FROZEN
			double AKE = 0;
			if ((SH2O + 0.0005) < SMC)
				AKE = SATRATIO;
			// UNFROZEN
			// RANGE OF VALIDITY FOR THE KERSTEN NUMBER (AKE)
			else {
				// KERSTEN NUMBER (USING "FINE" FORMULA, VALID FOR SOILS CONTAINING AT
				// LEAST 5% OF PARTICLES WITH DIAMETER LESS THAN 2.E-6 METERS.)
				// (FOR "COARSE" FORMULA, SEE PETERS-LIDARD ET AL., 1998).

				if (SATRATIO > 0.1) {

					AKE = Math.Log10(SATRATIO) + 1.0;
				}
				// USE K = KDRY
				else {

					AKE = 0.0;
				}
				//  THERMAL CONDUCTIVITY

			}

			return AKE * (THKSAT - THKDRY) + THKDRY;
//			if (double.IsNaN(DF))
//				throw new Exception("");
		}

		
		
		/// <summary>
		
		/// <summary>
		/// Respiration as a function of temperature
		/// </summary>
		/// <param name="T"></param>
		/// <returns></returns>
		public static double R(double T)
		{
//			if (T > 430)    //此句为刘永和自己加
//				T = 430;
			return Math.Exp(0.08 * (T - 298.16));
		}
		/// <summary>
		/// calculate soil water diffusivity and soil hydraulic conductivity.
		/// </summary>
		/// <param name="WDF"></param>
		/// <param name="WCND"></param>
		/// <param name="SMC"></param>
		/// <param name="FCR"></param>
		public static void WDFCND1(GridCell cell, out double WDF, out double WCND, double SMC, double FCR)
		{
			// soil water diffusivity
			double FACTR = Math.Max(0.01, SMC / cell.SMCMAX);
			double EXPON = cell.BEXP + 2.0;
			WDF = cell.DWSAT * Math.Pow(FACTR, EXPON);
			WDF *= 1.0 - FCR;

			// hydraulic conductivity
			EXPON = 2.0 * cell.BEXP + 3.0;
			WCND = cell.DKSAT * Math.Pow(FACTR, EXPON);
			WCND *= (1.0 - FCR);
		}
		
		/// <summary>
		/// calculate soil water diffusivity and soil hydraulic conductivity.
		/// 这个函数里的参数是造成地下水产流差别巨大的关键
		/// </summary>
		/// <param name="WDF">soil water diffusivity</param>
		/// <param name="WCND">hydraulic conductivity</param>
		/// <param name="SMC"></param>
		/// <param name="SICE"></param>
		public static void WDFCND2(GridCell cell, out double WDF, out double WCND, double SMC, double SICE)
		{
			// soil water diffusivity
			//当SMC(即外面的SH2O)的值大于0.3 m3/m3时，可能已经错了。后面WCND有问题的原因就在于这里
			double	FACTR = Math.Max(0.01, SMC / cell.SMCMAX);
			double	EXPON = cell.BEXP + 2.0;
			WDF = cell.DWSAT * Math.Pow(FACTR, EXPON);
			
			if (SICE > 0.0) {
				double	VKWGT = 1 / (1 + Math.Pow(500 * SICE, 3));
				WDF = VKWGT * WDF + (1 - VKWGT) * cell.DWSAT * Math.Pow(0.2 / cell.SMCMAX, EXPON);
			}
			// hydraulic conductivity
			EXPON = 2.0 * cell.BEXP + 3.0;
			WCND = cell.DKSAT * Math.Pow(FACTR, EXPON);
		}
		public static void Albedo(int VEGTYP, int  IST, int ISC, int ICE, int NSOIL,
			double DT, double COSZ, double ELAI, double ESAI,
			double TG, double TV, double SNOWH, double FSNO, double FWET,
			FortDoubleArray SMC, double SNEQVO, double SNEQV, double QSNOW, double FVEG,
			int ILOC, int JLOC,
			ref double ALBOLD, ref double TAUSS,//  //inout
			double[] ALBGRD, double[] ALBGRI, double[] ALBD, double[] ALBI, double[] FABD,//  //out
			double[] FABI, double[] FTDD, double[] FTID, double[] FTII, out double FSUN,//  //out
			double[] FREVI, double[] FREVD, double[] FREGD, double[] FREGI, out double BGAP,    //  //out
			out double WGAP)
		{


			// input
			//  INTEGER,                  INTENT(IN)  :: ILOC
			//  INTEGER,                  INTENT(IN)  :: JLOC
			//  INTEGER,                  INTENT(IN)  :: NSOIL  //number of soil layers
			//  INTEGER,                  INTENT(IN)  :: VEGTYP //vegetation type
			//  INTEGER,                  INTENT(IN)  :: IST    //surface type
			//  INTEGER,                  INTENT(IN)  :: ISC    //soil color type (1-lighest; 8-darkest)
			//  INTEGER,                  INTENT(IN)  :: ICE    //ice (ice = 1)
//
			//  REAL,                     INTENT(IN)  :: DT     //time step [sec]
			//  REAL,                     INTENT(IN)  :: QSNOW  //snowfall
			//  REAL,                     INTENT(IN)  :: COSZ   //cosine solar zenith angle for next time step
			//  REAL,                     INTENT(IN)  :: SNOWH  //snow height (mm)
			//  REAL,                     INTENT(IN)  :: TG     //ground temperature (k)
			//  REAL,                     INTENT(IN)  :: TV     //vegetation temperature (k)
			//  REAL,                     INTENT(IN)  :: ELAI   //LAI, one-sided, adjusted for burying by snow
			//  REAL,                     INTENT(IN)  :: ESAI   //SAI, one-sided, adjusted for burying by snow
			//  REAL,                     INTENT(IN)  :: FSNO   //fraction of grid covered by snow
			//  REAL,                     INTENT(IN)  :: FWET   //fraction of canopy that is wet
			//  REAL,                     INTENT(IN)  :: SNEQVO //snow mass at last time step(mm)
			//  REAL,                     INTENT(IN)  :: SNEQV  //snow mass (mm)
			//  REAL,                     INTENT(IN)  :: FVEG   //green vegetation fraction [0.0-1.0]
			//  REAL, DIMENSION(1:NSOIL), INTENT(IN)  :: SMC    //volumetric soil water (m3/m3)
//
			// inout
			//  REAL,                  INTENT(INOUT)  :: ALBOLD //snow albedo at last time step (CLASS type)
			//  REAL,                  INTENT(INOUT)  :: TAUSS  //non-dimensional snow age
//
			// output
			//  REAL, DIMENSION(1:    2), INTENT(OUT) :: ALBGRD //ground albedo (direct)
			//  REAL, DIMENSION(1:    2), INTENT(OUT) :: ALBGRI //ground albedo (diffuse)
			//  REAL, DIMENSION(1:    2), INTENT(OUT) :: ALBD   //surface albedo (direct)
			//  REAL, DIMENSION(1:    2), INTENT(OUT) :: ALBI   //surface albedo (diffuse)
			//  REAL, DIMENSION(1:    2), INTENT(OUT) :: FABD   //flux abs by veg (per unit direct flux)
			//  REAL, DIMENSION(1:    2), INTENT(OUT) :: FABI   //flux abs by veg (per unit diffuse flux)
			//  REAL, DIMENSION(1:    2), INTENT(OUT) :: FTDD   //down direct flux below veg (per unit dir flux)
			//  REAL, DIMENSION(1:    2), INTENT(OUT) :: FTID   //down diffuse flux below veg (per unit dir flux)
			//  REAL, DIMENSION(1:    2), INTENT(OUT) :: FTII   //down diffuse flux below veg (per unit dif flux)
			//  REAL,                     INTENT(OUT) :: FSUN   //sunlit fraction of canopy (-)
			//jref:start
			//  REAL, DIMENSION(1:    2), INTENT(OUT) :: FREVD
			//  REAL, DIMENSION(1:    2), INTENT(OUT) :: FREVI
			//  REAL, DIMENSION(1:    2), INTENT(OUT) :: FREGD
			//  REAL, DIMENSION(1:    2), INTENT(OUT) :: FREGI
			//  REAL, INTENT(OUT) :: BGAP
			//  REAL, INTENT(OUT) :: WGAP
			//  REAL                 :: FAGE     //snow age function
			double ALB = 0;
			//  INTEGER              :: IB       //indices
			//  INTEGER              :: NBAND    //number of solar radiation wave bands
			//  INTEGER              :: IC       //direct beam: ic=0; diffuse: ic=1
//
			//  REAL                 :: WL       //fraction of LAI+SAI that is LAI
			//  REAL                 :: WS       //fraction of LAI+SAI that is SAI
			//  REAL                 :: MPE      //prevents overflow for division by zero
//
			//  REAL, DIMENSION(1:2) :: RHO      //leaf/stem reflectance weighted by fraction LAI and SAI
			//  REAL, DIMENSION(1:2) :: TAU      //leaf/stem transmittance weighted by fraction LAI and SAI
			double[] FTDI = new double[2];     //down direct flux below veg per unit dif flux = 0
			double[] ALBSND = new double[2]; //   //snow albedo (direct)
			double[] ALBSNI = new double[2]; //   //snow albedo (diffuse)
//
			//  REAL                 :: VAI      //ELAI+ESAI
			double GDIR = 0;   //     //average projected leaf/stem area in solar direction
			//  REAL                 :: EXT      //optical depth direct beam per unit leaf + stem area
			// --------------------------------------------------------------------------------------------------

			int NBAND = 2;
			double MPE = 1E-06;
			BGAP = 0;
			WGAP = 0;
			FSUN = 0;
			//COSZ = 0;
			double[] RHO = new double[2];
			double[] TAU = new double[2];
			double[,] RHOL = Driver.vegparams.RHOL;  //leaf reflectance: 1=vis, 2=nir
			double[,] RHOS = Driver.vegparams.RHOS; //stem reflectance: 1=vis, 2=nir
			double[,] TAUL = Driver.vegparams.TAUL;// new double[30, 2]; //leaf transmittance: 1=vis, 2=nir
			double[,] TAUS = Driver.vegparams.TAUS; //stem transmittance: 1=vis, 2=nir
			double VAI = 0;      //ELAI+ESAI
			double WL = 0; //fraction of LAI+SAI that is LAI
			double WS = 0; //fraction of LAI+SAI that is SAI

			for (int IB = 0; IB < NBAND; IB++) {
				ALBD[IB] = 0;
				ALBI[IB] = 0;
				ALBGRD[IB] = 0;
				ALBGRI[IB] = 0;
				FABD[IB] = 0;
				FABI[IB] = 0;
				FTDD[IB] = 0;
				FTID[IB] = 0;
				FTII[IB] = 0;
				if (IB == 0)
					FSUN = 0;
			}

			if (COSZ <= 0)
				return;

			// weight reflectance/transmittance by LAI and SAI
			for (int IB = 0; IB < NBAND; IB++) {
				VAI = ELAI + ESAI;
				WL = ELAI / Math.Max(VAI, MPE);
				WS = ESAI / Math.Max(VAI, MPE);
				RHO[IB] = Math.Max(RHOL[VEGTYP - 1, IB] * WL + RHOS[VEGTYP - 1, IB] * WS, MPE);
				TAU[IB] = Math.Max(TAUL[VEGTYP - 1, IB] * WL + TAUS[VEGTYP - 1, IB] * WS, MPE);
			}
			if (double.IsNaN(VAI))
				throw new Exception();
			// snow age
			double FAGE = 0;
			SNOW_AGE(DT, TG, SNEQVO, SNEQV, ref TAUSS, out FAGE);

			//snow albedos: only if COSZ > 0 and FSNO > 0
			if (NoahMP.OPT_ALB == 1) {
				SNOWALB_BATS(NBAND, FSNO, COSZ, FAGE, ALBSND, ALBSNI);
			}
			if (NoahMP.OPT_ALB == 2) {
				SNOWALB_CLASS(NBAND, QSNOW, DT, ref ALB, ALBOLD, ALBSND, ALBSNI);
				ALBOLD = ALB;
			}

			// ground surface albedo
			GROUNDALB(NSOIL, NBAND, ICE, IST, ISC,
				FSNO, SMC, ALBSND, ALBSNI, COSZ,
				TG, ILOC, JLOC,
				ALBGRD, ALBGRI);                              //out

			// loop over NBAND wavebands to calculate surface albedos and solar
			// fluxes for unit incoming direct (IC=0) and diffuse flux (IC=1)

			for (int IB = 0; IB < NBAND; IB++) {
				int IC = 0;      // direct

				
				//在这个函数中，ALBGRD只是输出变量，不修改
				TWOSTREAM(IB, IC, VEGTYP, COSZ, VAI,
					FWET, TV, ALBGRD, ALBGRI, RHO,
					TAU, FVEG, IST, ILOC, JLOC,
					FABD, ALBD, FTDD, FTID, out GDIR,//)   //out
					FREVD, FREGD, out BGAP, out WGAP);

				IC = 1;      // diffuse
				TWOSTREAM(IB, IC, VEGTYP, COSZ, VAI,
					FWET, TV, ALBGRD, ALBGRI, RHO,
					TAU, FVEG, IST, ILOC, JLOC,
					FABI, ALBI, FTDI, FTII, out GDIR, //   //out
					FREVI, FREGI, out BGAP, out WGAP);

			}

			// sunlit fraction of canopy. set FSUN = 0 if FSUN < 0.01.
			double EXT = GDIR / COSZ * Math.Sqrt(1 - RHO[0] - TAU[0]);
			FSUN = (1 - Math.Exp(-EXT * VAI)) / Math.Max(EXT * VAI, MPE);
			EXT = FSUN;

			if (EXT < 0.01) {
				WL = 0;
			} else {
				WL = EXT;
			}
			FSUN = WL;


		}
		public static void SNOW_AGE(double DT, double TG, double SNEQVO, double SNEQV, ref double  TAUSS, out double FAGE)
		{

			//input
			//   REAL, INTENT(IN) :: DT        //main time step (s)
			//   REAL, INTENT(IN) :: TG        //ground temperature (k)
			//   REAL, INTENT(IN) :: SNEQVO    //snow mass at last time step(mm)
			//   REAL, INTENT(IN) :: SNEQV     //snow water per unit ground area (mm)
//
			//output
			//   REAL, INTENT(OUT) :: FAGE     //snow age
//
			//input/output
			//   REAL, INTENT(INOUT) :: TAUSS      //non-dimensional snow age
			//local
			double TAGE;       //total aging effects
			double AGE1;      //effects of grain growth due to vapor diffusion
			double AGE2;      //effects of grain growth at freezing of melt water
			double AGE3;      //effects of soot
			double DELA;      //temporary variable
			double SGE;      //temporary variable
			double DELS;      //temporary variable
			double DELA0;    //temporary variable
			double ARG;       //temporary variable
			// See Yang et al. (1997) J.of Climate for detail.
			//---------------------------------------------------------------------------------------------------
			
			if (SNEQV <= 0.0) {
				TAUSS = 0;
			} else if (SNEQV > 800) {
				TAUSS = 0;
			} else {
				DELA0 = 1E-6 * DT;
				ARG = 5E3 * (1 / NoahMP.TFRZ - 1 / TG);
				AGE1 = Math.Exp(ARG);
				AGE2 = Math.Exp(Math.Min(0, 10 * ARG));
				AGE3 = 0.3;
				TAGE = AGE1 + AGE2 + AGE3;
				DELA = DELA0 * TAGE;
				DELS = Math.Max(0.0, SNEQV - SNEQVO) / NoahMP.SWEMX;
				SGE = (TAUSS + DELA) * (1.0 - DELS);
				TAUSS = Math.Max(0, SGE);
			}

			FAGE = TAUSS / (TAUSS + 1);

		}
		public static void GROUNDALB(int NSOIL, int NBAND, int ICE, int IST, int ISC,
			double FSNO, FortDoubleArray SMC, double[] ALBSND, double[] ALBSNI, double COSZ,
			double  TG, int ILOC, int JLOC,
			double[]  ALBGRD, double[] ALBGRI)                              //out
		{
			//input

			//  INTEGER,                  INTENT(IN)  :: ILOC   //grid index
			//  INTEGER,                  INTENT(IN)  :: JLOC   //grid index
			//  INTEGER,                  INTENT(IN)  :: NSOIL  //number of soil layers
			//  INTEGER,                  INTENT(IN)  :: NBAND  //number of solar radiation waveband classes
			//  INTEGER,                  INTENT(IN)  :: ICE    //value of ist for land ice
			//  INTEGER,                  INTENT(IN)  :: IST    //surface type
			//  INTEGER,                  INTENT(IN)  :: ISC    //soil color class (1-lighest; 8-darkest)
			//  REAL,                     INTENT(IN)  :: FSNO   //fraction of surface covered with snow (-)
			//  REAL,                     INTENT(IN)  :: TG     //ground temperature (k)
			//  REAL,                     INTENT(IN)  :: COSZ   //cosine solar zenith angle (0-1)
			//  REAL, DIMENSION(1:NSOIL), INTENT(IN)  :: SMC    //volumetric soil water content (m3/m3)
			//  REAL, DIMENSION(1:    2), INTENT(IN)  :: ALBSND //direct beam snow albedo (vis, nir)
			//  REAL, DIMENSION(1:    2), INTENT(IN)  :: ALBSNI //diffuse snow albedo (vis, nir)

			//output

			//double[] ALBGRD;  //ground albedo (direct beam: vis, nir)
			//double[] ALBGRI;  //ground albedo (diffuse: vis, nir)

			//local

			//int IB;//     //waveband number (1=vis, 2=nir)
			double INC;//    //soil water correction factor for soil albedo
			double ALBSOD; //soil albedo (direct)
			double ALBSOI; //soil albedo (diffuse)
			// --------------------------------------------------------------------------------------------------

			for (int IB = 0; IB < NBAND; IB++) {

				INC = Math.Max(0.11 - 0.40 * SMC[1], 0);
				if (IST == 1) {                       //soil
					ALBSOD = Math.Min(RAD_PARAMS.ALBSAT[ISC - 1, IB] + INC, RAD_PARAMS.ALBDRY[ISC - 1, IB]);
					ALBSOI = ALBSOD;
				} else if (TG > NoahMP.TFRZ) {                //unfrozen lake, wetland
					ALBSOD = 0.06 / (Math.Pow(Math.Max(0.01, COSZ), 1.7) + 0.15);
					ALBSOI = 0.06;
				} else {                                      //frozen lake, wetland
					ALBSOD = RAD_PARAMS.ALBLAK[IB];
					ALBSOI = ALBSOD;
				}

				//increase desert and semi-desert albedos

				if (IST == 1 && ISC == 9) {
					ALBSOD += 0.10;
					ALBSOI += 0.10;
				}

				ALBGRD[IB] = ALBSOD * (1 - FSNO) + ALBSND[IB] * FSNO;
				ALBGRI[IB] = ALBSOI * (1 - FSNO) + ALBSNI[IB] * FSNO;
				
			}
		}
		
		public static void SNOWALB_BATS(int NBAND, double FSNO, double COSZ, double FAGE, double[] ALBSND, double[] ALBSNI)
		{
			// input
//
			//  INTEGER,INTENT(IN) :: NBAND  //number of waveband classes
//
			//  REAL,INTENT(IN) :: COSZ    //cosine solar zenith angle
			//  REAL,INTENT(IN) :: FSNO    //snow cover fraction (-)
			//  REAL,INTENT(IN) :: FAGE    //snow age correction
//
			// output
//
			//  REAL, DIMENSION(1:2),INTENT(OUT) :: ALBSND //snow albedo for direct(1=vis, 2=nir)
			//  REAL, DIMENSION(1:2),INTENT(OUT) :: ALBSNI //snow albedo for diffuse
			// ---------------------------------------------------------------------------------------------
//
			// ------------------------ local variables ----------------------------------------------------
			//  INTEGER :: IB          //waveband class
//
			double FZEN;                 //zenith angle correction
			double CF1;                 //temperary variable
			double SL2;                //2.0*SL
			double SL1;                //1/SL
			double SL;  //                   //adjustable parameter
			double C1 = 0.2; //  //default in BATS
			double C2 = 0.5; //  //default in BATS
			//double C1 = 0.2 * 2;// // double the default to match Sleepers River's
			//double C2 = 0.5 * 2; //. // snow surface albedo (double aging effects)
			// ---------------------------------------------------------------------------------------------
			// zero albedos for all points

			for (int i = 0; i < NBAND; i++) {
				ALBSND[i] = 0;
				ALBSNI[i] = 0;
			}
			// when cosz > 0

			SL = 2.0;
			SL1 = 1 / SL;
			SL2 = 2 * SL;
			CF1 = ((1 + SL1) / (1 + SL2 * COSZ) - SL1);
			FZEN = Math.Max(CF1, 0);

			ALBSNI[0] = 0.95 * (1 - C1 * FAGE);
			ALBSNI[1] = 0.65 * (1 - C2 * FAGE);
			ALBSND[0] = ALBSNI[0] + 0.4 * FZEN * (1 - ALBSNI[0]);    //  vis direct
			ALBSND[1] = ALBSNI[1] + 0.4 * FZEN * (1 - ALBSNI[1]);    // nir direct
			
		}

		public static void SNOWALB_CLASS(int NBAND, double QSNOW, double DT, ref double ALB, double ALBOLD, double[]  ALBSND, double[] ALBSNI)
		{
			// input

			//  INTEGER,INTENT(IN) :: ILOC //grid index
			//  INTEGER,INTENT(IN) :: JLOC //grid index
			//  INTEGER,INTENT(IN) :: NBAND  //number of waveband classes
//
			//  REAL,INTENT(IN) :: QSNOW     //snowfall (mm/s)
			//  REAL,INTENT(IN) :: DT        //time step (sec)
			//  REAL,INTENT(IN) :: ALBOLD    //snow albedo at last time step

			// in  out

			//  REAL,                INTENT(INOUT) :: ALB        // 
			// output

			// REAL, DIMENSION(1:2),INTENT(OUT) :: ALBSND //snow albedo for direct(1=vis, 2=nir)
			// REAL, DIMENSION(1:2),INTENT(OUT) :: ALBSNI //snow albedo for diffuse
			// ---------------------------------------------------------------------------------------------

			// ------------------------ local variables ----------------------------------------------------
			int IB;          //waveband class

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

			if (QSNOW > 0) {
				ALB = ALB + Math.Min(QSNOW * DT, NoahMP.SWEMX) * (0.84 - ALB) / (NoahMP.SWEMX);
			}
//			if (ALB == 0)
//				throw new Exception("");

			ALBSNI[0] = ALB;         // vis diffuse
			ALBSNI[1] = ALB;        // nir diffuse
			ALBSND[0] = ALB;        // vis direct
			ALBSND[1] = ALB;        // nir direct
		}
		public static	double F1_Func(double AB, double BC)
		{
			return Math.Pow(AB, (BC - 25) / 10);
		}
		public static	double F2_Func(double AB)
		{
			return 1 + Math.Exp((-2.2E05 + 710 * (AB + 273.16)) / (8.314 * (AB + 273.16)));
		}
		

	}
}
