/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2018/1/20
 * Time: 16:18
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;

namespace NoahMP
{
	/// <summary>
	/// Description of SFCDIF1.
	/// </summary>
	public static class SFCDIF1
	{
	  public static	void SFCDIF(int ITER, double SFCTMP, double RHOAIR, double H, double QAIR,  //in
			double ZLVL, double ZPD, double Z0M, double Z0H, double UR,  //in
			double  MPE, int ILOC, int JLOC,                  //in
			ref  double  MOZ, ref double MOZSGN, ref double FM, ref double FH, ref double FM2, ref double FH2,  //inout
			out double CM, out double  CH, ref double FV, out double CH2) //out
		{
			// -------------------------------------------------------------------------------------------------
			// computing surface drag coefficient CM for momentum and CH for heat
			// -------------------------------------------------------------------------------------------------
			//   IMPLICIT NONE
			// -------------------------------------------------------------------------------------------------
			// inputs

			//    INTEGER,              INTENT(IN) :: ILOC   //grid index
			//    INTEGER,              INTENT(IN) :: JLOC   //grid index
			//    INTEGER,              INTENT(IN) :: ITER   //iteration index
			//    REAL,                 INTENT(IN) :: SFCTMP //temperature at reference height (k)
			//    REAL,                 INTENT(IN) :: RHOAIR //density air (kg/m**3)
			//    REAL,                 INTENT(IN) :: H      //sensible heat flux (w/m2) [+ to atm]
			//    REAL,                 INTENT(IN) :: QAIR   //specific humidity at reference height (kg/kg)
			//    REAL,                 INTENT(IN) :: ZLVL   //reference height  (m)
			//    REAL,                 INTENT(IN) :: ZPD    //zero plane displacement (m)
			//    REAL,                 INTENT(IN) :: Z0H    //roughness length, sensible heat, ground (m)
			//    REAL,                 INTENT(IN) :: Z0M    //roughness length, momentum, ground (m)
			//    REAL,                 INTENT(IN) :: UR     //wind speed (m/s)
			//    REAL,                 INTENT(IN) :: MPE    //prevents overflow error if division by zero
			//// in  out
//
			//    INTEGER,           INTENT(INOUT) :: MOZSGN //number of times moz changes sign
			//    REAL,              INTENT(INOUT) :: MOZ    //Monin-Obukhov stability (z/L)
			//    REAL,              INTENT(INOUT) :: FM     //momentum stability correction, weighted by prior iters
			//    REAL,              INTENT(INOUT) :: FH     //sen heat stability correction, weighted by prior iters
//		 REAL,              INTENT(INOUT) :: FM2    //sen heat stability correction, weighted by prior iters
			//    REAL,              INTENT(INOUT) :: FH2    //sen heat stability correction, weighted by prior iters
//
			//// outputs
//
			//    REAL,                INTENT(OUT) :: CM     //drag coefficient for momentum
			//    REAL,                INTENT(OUT) :: CH     //drag coefficient for heat
			//    REAL,                INTENT(OUT) :: FV     //friction velocity (m/s)
			//    REAL,                INTENT(OUT) :: CH2    //drag coefficient for heat
//
			//// locals
			double MOL;                      //Monin-Obukhov length (m)
			//    REAL    :: TMPCM                    //temporary calculation for CM
			//    REAL    :: TMPCH                    //temporary calculation for CH
			double FMNEW;                    //stability correction factor, momentum, for current moz
			double FHNEW;                    //stability correction factor, sen heat, for current moz
			//    REAL    :: MOZOLD                   //Monin-Obukhov stability parameter from prior iteration
			//    REAL    :: TMP1,TMP2,TMP3,TMP4,TMP5 //temporary calculation
			double TVIR;                     //temporary virtual temperature (k)
			double MOZ2;                     //2/L
			//    REAL    :: TMPCM2                   //temporary calculation for CM2
			//    REAL    :: TMPCH2                   //temporary calculation for CH2
			double FM2NEW;                   //stability correction factor, momentum, for current moz
			double FH2NEW;                   //stability correction factor, sen heat, for current moz
			double TMP12, TMP22, TMP32;       //temporary calculation

			double CMFM, CHFH, CM2FM2, CH2FH2;
			// -------------------------------------------------------------------------------------------------
			// Monin-Obukhov stability parameter moz for next iteration

			double MOZOLD = MOZ;

			if (ZLVL <= ZPD) {
				//write(*,*) 'critical problem: ZLVL <= ZPD; model stops'
				//call wrf_error_fatal("STOP in Noah-MP")
			}

			double	TMPCM = Math.Log((ZLVL - ZPD) / Z0M);
			double	TMPCH = Math.Log((ZLVL - ZPD) / Z0H);
			double TMPCM2 = Math.Log((2.0 + Z0M) / Z0M);
			double TMPCH2 = Math.Log((2.0 + Z0H) / Z0H);
			
			if (ITER == 1) {
				FV = 0.0;
				MOZ = 0.0;
				MOL = 0.0;
				MOZ2 = 0.0;
			} else {
				TVIR = (1 + 0.61 * QAIR) * SFCTMP;
				double TMP1 = NoahMP.VKC * (NoahMP.GRAV / TVIR) * H / (RHOAIR * NoahMP.CPAIR);
				if (Math.Abs(TMP1) <= MPE)
					TMP1 = MPE;
				MOL = -1 * Math.Pow(FV, 3) / TMP1;
				MOZ = Math.Min((ZLVL - ZPD) / MOL, 1);
				MOZ2 = Math.Min((2.0 + Z0H) / MOL, 1);
//				if (double.IsNaN(MOZ2))
//					throw new Exception();
			}
            
			// accumulate number of times moz changes sign.

			if (MOZOLD * MOZ < 0) {
				MOZSGN = MOZSGN + 1;
			}
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
				double TMP1 = Math.Pow(1 - 16 * MOZ, 0.25);
				double TMP2 = Math.Log((1 + TMP1 * TMP1) / 2);
				double TMP3 = Math.Log((1 + TMP1) / 2);
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
			double VKC = NoahMP.VKC;
			CM = VKC * VKC / (CMFM * CMFM);
			CH = VKC * VKC / (CMFM * CHFH);
			CH2 = VKC * VKC / (CM2FM2 * CH2FH2);
			// friction velocity

			FV = UR * Math.Sqrt(CM);
			
			CH2 = VKC * FV / CH2FH2;
//			if (double.IsNaN(CH) || double.IsNaN(FV))
//				throw new Exception();

		}
	}
}
