/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2018/1/10
 * Time: 10:27
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;

namespace NoahMP
{
	public class WATER_Input
	{
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
	}
}


