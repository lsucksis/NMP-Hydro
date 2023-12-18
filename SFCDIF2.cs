/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2018/1/20
 * Time: 14:17
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;

namespace NoahMP
{
	/// <summary>
	/// Description of SFCDIF.
	/// </summary>
	public static class SFCDIF2
	{
		// LECH'S SURFACE FUNCTIONS
		public static double PSLMU(double ZZ)
		{
			return -0.96 * Math.Log(1.0 - 4.5 * ZZ);
		}
		public static double PSLMS(double ZZ)
		{
			return ZZ * RRIC - 2.076 * (1 - 1.0 / (ZZ + 1));
		}
		public static double PSLHU(double ZZ)
		{
			return -0.96 * Math.Log(1.0 - 4.5 * ZZ);
		}
		public static double PSLHS(double ZZ)
		{
			return ZZ * RFAC - 2.076 * (1 - 1 / (ZZ + 1));
		}
		// PAULSON'S SURFACE FUNCTIONS
		public static double PSPMU(double XX)
		{
			return -2 * Math.Log((XX + 1) * 0.5) - Math.Log((XX * XX + 1) * 0.5) + 2 * Math.Atan(XX) - PIHF;
		}
		public static double PSPMS(double YY)
		{
			return 5 * YY;
		}
		public static double PSPHU(double XX)
		{
			return -2 * Math.Log((XX * XX + 1) * 0.5);
		}
		public static double PSPHS(double YY)
		{
			return 5 * YY;
		}
		static int ITRMX = 5;
		static double WWST = 1.2;
		static double WWST2 = WWST * WWST;
		static double VKRM = 0.40;
		static double EXCM = 0.001;
		static double BETA = 1.0 / 270.0;
		static double BTG = BETA * NoahMP.GRAV;
		static double ELFC = VKRM * BTG;
		static double WOLD = 0.15;
		static double WNEW = 1.0 - WOLD;
		static double PIHF = 3.14159265 / 2;
		static double EPSU2 = 1e-4;
		static double EPSUST = 0.07;
		static double EPSIT = 1e-4;
		static double EPSA = 1e-8;
		static double ZTMIN = -5.0;
		static double ZTMAX = 1.0;
		static double HPBL = 1000.0;
		static double SQVISC = 258.2;
		static double RIC = 0.183;
		static double RRIC = 1.0 / RIC;
		static double FHNEU = 0.8;
		static double RFC = 0.191;
		static double RFAC = RIC / (FHNEU * RFC * RFC);
		public static void SFCDIF(int ITER, double  Z0, double THZ0, double THLM, double SFCSPD,  //in
			double CZIL, double ZLM, double ILOC, double JLOC,          //in
			ref  double AKMS, ref  double AKHS, ref  double RLMO, ref double WSTAR2,         //inout
			ref  double USTAR)                                   //inout  原有Noah-MP代码可能不对
		{

// -------------------------------------------------------------------------------------------------
// SUBROUTINE SFCDIF (renamed SFCDIF_off to avoid clash with Eta PBL)
// -------------------------------------------------------------------------------------------------
// CALCULATE SURFACE LAYER EXCHANGE COEFFICIENTS VIA ITERATIVE PROCESS.
// SEE CHEN ET AL (1997, BLM)
// -------------------------------------------------------------------------------------------------

//    INTEGER, INTENT(IN) :: ILOC
//    INTEGER, INTENT(IN) :: JLOC
//    INTEGER, INTENT(IN) :: ITER
//    REAL,    INTENT(IN) :: ZLM, Z0, THZ0, THLM, SFCSPD, CZIL
//    REAL, intent(INOUT) :: AKMS
//    REAL, intent(INOUT) :: AKHS
//    REAL, intent(INOUT) :: RLMO
//    REAL, intent(INOUT) :: WSTAR2
//    REAL,   intent(OUT) :: USTAR
//
//    REAL     ZZ, PSLMU, PSLMS, PSLHU, PSLHS
//    REAL     XX, PSPMU, YY, PSPMS, PSPHU, PSPHS
//    REAL     ZILFC, ZU, ZT, RDZ, CXCH
//    REAL     DTHV, DU2, BTGH, ZSLU, ZSLT, RLOGU, RLOGT
//    REAL     ZETALT, ZETALU, ZETAU, ZETAT, XLU4, XLT4, XU4, XT4
//
			double PSMZ, SIMM, PSHZ, SIMH, USTARK, RLMN, RLMA;

//    INTEGER  ILECH, ITR
		
			// ----------------------------------------------------------------------
// NOTE: THE TWO CODE BLOCKS BELOW DEFINE FUNCTIONS
// ----------------------------------------------------------------------


// THIS ROUTINE SFCDIF CAN HANDLE BOTH OVER OPEN WATER (SEA, OCEAN) AND
// OVER SOLID SURFACE (LAND, SEA-ICE).
// ----------------------------------------------------------------------
//     ZTFC: RATIO OF ZOH/ZOM  LESS OR EQUAL THAN 1
//     C......ZTFC=0.1
//     CZIL: CONSTANT C IN Zilitinkevich, S. S.1995,:NOTE ABOUT ZT
// ----------------------------------------------------------------------
			int ILECH = 0;

// ----------------------------------------------------------------------
			double ZILFC = -CZIL * VKRM * SQVISC;
			double ZU = Z0;
			double RDZ = 1 / ZLM;
			double CXCH = EXCM * RDZ;
			double DTHV = THLM - THZ0;

// BELJARS CORRECTION OF USTAR
			double DU2 = Math.Max(SFCSPD * SFCSPD, EPSU2);
			double BTGH = BTG * HPBL;
			
			//USTAR=0;
			if (ITER == 1) {
				if (BTGH * AKHS * DTHV != 0.0) {
					double temp = Math.Abs(BTGH * AKHS * DTHV);
					WSTAR2 = WWST2 * Math.Pow(temp, 2.0 / 3);
				} else {
					WSTAR2 = 0.0;
				}
				USTAR = Math.Max(Math.Sqrt(AKMS * Math.Sqrt(DU2 + WSTAR2)), EPSUST);
				RLMO = ELFC * AKHS * DTHV / Math.Pow(USTAR, 3);
			}

// ZILITINKEVITCH APPROACH FOR ZT
			double ZT = Math.Max(1e-6, Math.Exp(ZILFC * Math.Sqrt(USTAR * Z0)) * Z0);
			double ZSLU = ZLM + ZU;
			double ZSLT = ZLM + ZT;
			double RLOGU = Math.Log(ZSLU / ZU);
			double RLOGT = Math.Log(ZSLT / ZT);
// 1./MONIN-OBUKKHOV LENGTH-SCALE
// ----------------------------------------------------------------------
			double ZETALT = Math.Max(ZSLT * RLMO, ZTMIN);
			RLMO = ZETALT / ZSLT;
			double ZETALU = ZSLU * RLMO;
			double ZETAU = ZU * RLMO;
			double ZETAT = ZT * RLMO;
    	
			if (ILECH == 0) {
				if (RLMO < 0) {
					double XLU4 = 1 - 16 * ZETALU;
					double	XLT4 = 1 - 16.0 * ZETALT;
					double XU4 = 1 - 16.0 * ZETAU;
					double	XT4 = 1 - 16.0 * ZETAT;
					double	XLU = Math.Sqrt(Math.Sqrt(XLU4));
					double	XLT = Math.Sqrt(Math.Sqrt(XLT4));
					double	XU = Math.Sqrt(Math.Sqrt(XU4));

					double	XT = Math.Sqrt(Math.Sqrt(XT4));
					PSMZ = PSPMU(XU);
					SIMM = PSPMU(XLU) - PSMZ + RLOGU;
					PSHZ = PSPHU(XT);
					SIMH = PSPHU(XLT) - PSHZ + RLOGT;
				} else {
					ZETALU = Math.Min(ZETALU, ZTMAX);
					ZETALT = Math.Min(ZETALT, ZTMAX);
					PSMZ = PSPMS(ZETAU);
					SIMM = PSPMS(ZETALU) - PSMZ + RLOGU;
					PSHZ = PSPHS(ZETAT);
					SIMH = PSPHS(ZETALT) - PSHZ + RLOGT;
				}
// ----------------------------------------------------------------------
// LECH'S FUNCTIONS
// ----------------------------------------------------------------------
			} else {
				if (RLMO < 0) {
					PSMZ = PSLMU(ZETAU);
					SIMM = PSLMU(ZETALU) - PSMZ + RLOGU;
					PSHZ = PSLHU(ZETAT);
					SIMH = PSLHU(ZETALT) - PSHZ + RLOGT;
				} else {
					ZETALU = Math.Min(ZETALU, ZTMAX);
					ZETALT = Math.Min(ZETALT, ZTMAX);
					PSMZ = PSLMS(ZETAU);
					SIMM = PSLMS(ZETALU) - PSMZ + RLOGU;
					PSHZ = PSLHS(ZETAT);
					SIMH = PSLHS(ZETALT) - PSHZ + RLOGT;
				}
// ----------------------------------------------------------------------
			}

// ----------------------------------------------------------------------
// BELJAARS CORRECTION FOR USTAR
// ----------------------------------------------------------------------
			USTAR = Math.Max(Math.Sqrt(AKMS * Math.Sqrt(DU2 + WSTAR2)), EPSUST);

// ZILITINKEVITCH FIX FOR ZT
			ZT = Math.Max(1e-6, Math.Exp(ZILFC * Math.Sqrt(USTAR * Z0)) * Z0);
			ZSLT = ZLM + ZT;
//-----------------------------------------------------------------------
			RLOGT = Math.Log(ZSLT / ZT);
			USTARK = USTAR * VKRM;
			AKMS = Math.Max(USTARK / SIMM, CXCH);
//-----------------------------------------------------------------------
// if STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO
//-----------------------------------------------------------------------
			AKHS = Math.Max(USTARK / SIMH, CXCH);

			if (BTGH * AKHS * DTHV != 0.0) {
				WSTAR2 = WWST2 * Math.Pow(Math.Abs(BTGH * AKHS * DTHV), 2.0 / 3.0);
			} else {
				WSTAR2 = 0.0;
			}
//-----------------------------------------------------------------------
			RLMN = ELFC * AKHS * DTHV / Math.Pow(USTAR, 3);
//-----------------------------------------------------------------------
//     if(ABS((RLMN-RLMO)/RLMA)<EPSIT)    GO TO 110
//-----------------------------------------------------------------------
			RLMA = RLMO * WOLD + RLMN * WNEW;
//-----------------------------------------------------------------------
			RLMO = RLMA;

//       write(*,'(a20,10f15.6)')'SFCDif: RLMO=',RLMO,RLMN,ELFC , AKHS , DTHV , USTAR
//    END DO
// ----------------------------------------------------------------------
		}
		
	}
}
