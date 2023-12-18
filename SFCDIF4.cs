/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2018/1/20
 * Time: 16:20
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;

namespace NoahMP
{
	/// <summary>
	/// Description of SFCDIF4.
	/// </summary>
	public static class SFCDIF4
	{
		public static void   SFCDIF(int ILOC, int JLOC, double UX, double VX, double T1D,
			double P1D, double PSFCPA, double PBLH, double DX, double ZNT, double
		                      TSK, double QX, double ZLVL, double IZ0TLND, ref double QSFC, double
		                      HFX, double QFX, ref double CM, ref double CHS, ref double CHS2, ref double
		                      CQS2, ref  double RMOL, ref  double UST, out double U10, out double V10)
		{

			
			//  Compute surface drag coefficients CM for momentum and CH for heat
			//  Joakim Refslund, 2011. Modified from YSU SFCLAY.
			
			// parameters
			double XKA = 2.4E-5;
			double PRT = 1;     //prandtl number

			// input

			//   REAL,   INTENT(IN )   :: PBLH      // planetary boundary layer height
			//   REAL,   INTENT(IN )   :: TSK       // skin temperature
			//   REAL,   INTENT(IN )   :: PSFCPA    // pressure in pascal
			//   REAL,   INTENT(IN )   :: P1D       //lowest model layer pressure (Pa)
			//   REAL,   INTENT(IN )   :: T1D       //lowest model layer temperature
			//   REAL,   INTENT(IN )   :: QX        //water vapor mixing ratio (kg/kg)
			//   REAL,   INTENT(IN )   :: QX        //water vapor specific humidity (kg/kg)
			//   REAL,   INTENT(IN )   :: ZLVL      // thickness of lowest full level layer
			//   REAL,   INTENT(IN )   :: HFX       // sensible heat flux
			//   REAL,   INTENT(IN )   :: QFX       // moisture flux
			//   REAL,   INTENT(IN )   :: DX        // horisontal grid spacing
			//   REAL,   INTENT(IN )   :: UX
//			REAL,   INTENT(IN )   :: VX
			//   REAL,   INTENT(IN )   :: ZNT
			//   REAL,   INTENT(INOUT ) :: QSFC

			//   REAL,   INTENT(INOUT) :: RMOL
			//   REAL,   INTENT(INOUT) :: UST
			//   REAL,   INTENT(INOUT) :: CHS2
			//   REAL,   INTENT(INOUT) :: CQS2
			//   REAL,   INTENT(INOUT) :: CHS
			//   REAL,   INTENT(INOUT) :: CM

			// diagnostics out
			//   REAL,   INTENT(OUT)   :: U10
			//   REAL,   INTENT(OUT)   :: V10
			//   REAL,   INTENT(OUT)   :: TH2
			//   REAL,   INTENT(OUT)   :: T2
			//   REAL,   INTENT(OUT)   :: Q2
			//   REAL,   INTENT(OUT)   :: QSFC

			// optional vars
			//   INTEGER,OPTIONAL,INTENT(IN ) :: IZ0TLND
//			// local
			//   INTEGER :: REGIME  // Stability regime
			//   REAL    :: ZA      // Height of full-sigma level
			//   REAL    :: THVX    // Virtual potential temperature
			//   REAL    :: ZQKL    // Height of upper half level
			//   REAL    :: ZQKLP1  // Height of lower half level (surface)
			//   REAL    :: THX     // Potential temperature
			double PSIH = 0;   // similarity function for heat
			double PSIH2 = 0;  // Similarity function for heat 2m
			double PSIH10 = 0;  // Similarity function for heat 10m
			double PSIM = 0;   // similarity function for momentum
			//   REAL    :: PSIM2   // Similarity function for momentum 2m
			double PSIM10 = 0;  // Similarity function for momentum 10m
			//   REAL    :: DENOMQ  // Denominator used for flux calc.
			//   REAL    :: DENOMQ2 // Denominator used for flux calc.
			//   REAL    :: DENOMT2 // Denominator used for flux calc.
			//   REAL    :: WSPDI   // Initial wind speed
			//   REAL    :: GZ1OZ0  // log(za/z0)
			//   REAL    :: GZ2OZ0  // log(z2/z0)
			//   REAL    :: GZ10OZ0 // log(z10/z0)
			//   REAL    :: RHOX    // density
			//   REAL    :: GOVRTH  // g/theta for stability L
			//   REAL    :: TGDSA   // tsk
			//   REAL    :: SCR3    // temporal variable -> input variable T1D
			//   REAL    :: TVIR    // temporal variable SRC4 -> TVIR
			//   REAL    :: THGB    // Potential temperature ground
			//   REAL    :: PSFC    // Surface pressure
			//   REAL    :: BR      // bulk richardson number
			//   REAL    :: CPM
			//   REAL    :: MOL
			//   REAL    :: ZOL
			//   REAL    :: QGH
//			   REAL    :: WSPD
//
			//   INTEGER :: N,I,K,KK,L,NZOL,NK,NZOL2,NZOL10
//
			//   REAL    ::  PL,THCON,TVCON,E1
			//   REAL    ::  ZL,TSKV,DTHVDZ,DTHVM,VCONV,RZOL,RZOL2,RZOL10,ZOL2,ZOL10
			//   REAL    ::  DTG,PSIX,DTTHX,PSIX10,PSIT,PSIT2,PSIQ,PSIQ2,PSIQ10
			//   REAL    ::  FLUXC,VSGD,Z0Q,VISC,RESTAR,CZIL,RESTAR2
			//-------------------------------------------------------------------
			double PSIT = 0;
			double MOL = 1 / RMOL;
			double ZL = 0.01;
			double PSFC = PSFCPA / 1000;
			double P1000mb = 100000;

			// convert (tah or tgb = tsk) temperature to potential temperature.
			double TGDSA = TSK;
			double THGB = TSK * Math.Pow(P1000mb / PSFCPA, Constants.rcp);

			// store virtual, virtual potential and potential temperature
			double PL = P1D / 1000;
			double THX = T1D * Math.Pow(P1000mb * 0.001 / PL, Constants.rcp);
			double THVX = THX * (1 + Constants.EP_1 * QX);
			double TVIR = T1D * (1 + Constants.EP_1 * QX);

			// for land points QSFC can come from previous time step
			//QSFC=EP_2*E1/(PSFC-E1)
			
			if (QSFC <= 0.0) {
				//testing this
				double E1 = Constants.SVP1 * Math.Exp(Constants.SVP2 * (TGDSA - Constants.SVPT0) / (TGDSA - Constants.SVP3));
				QSFC = Constants.EP_2 * E1 / (PSFC - E1);
				if (double.IsNaN(QSFC))
					throw new Exception("");
				//write(*,*) "JREF: IN SFCDIF4, QSFC WAS NEG. NOW = ",QSFC
			}
			// qgh changed to use lowest-level air temp consistent with myjsfc change
			// q2sat = qgh in lsm
			//jref: canres and esat is calculated in the loop so should that be changed??
			//   QGH=EP_2*E1/(PL-E1)
			double CPM = Constants.cp * (1 + 0.8 * QX);

			// compute the height of half-sigma levels above ground level
			//ZA=0.5*DZ8W
			double ZA = ZLVL;

			// compute density and part of monin-obukhov length L
			double RHOX = PSFC * 1000 / (Constants.R_d * TVIR);
			double GOVRTH = Constants.g / THX;

			// calculate bulk richardson no. of surface layer,
			// according to akb(1976), eq(12).
			double GZ1OZ0 = Math.Log(ZA / ZNT);
			double GZ2OZ0 = Math.Log(2 / ZNT);
			double GZ10OZ0 = Math.Log(10 / ZNT);
			double WSPD = Math.Sqrt(UX * UX + VX * VX);

			// virtual pot. temperature difference between input layer and lowest model layer
			double TSKV = THGB * (1 + Constants.EP_1 * QSFC);
			double DTHVDZ = (THVX - TSKV);
			// convective velocity scale Vc and subgrid-scale velocity Vsg
			// following Beljaars (1995, QJRMS) and Mahrt and Sun (1995, MWR)
			//                                ... HONG Aug. 2001
//
			// VCONV = 0.25*Sqrt(g/tskv*pblh(i)*dthvm)
			// use Beljaars over land, old MM5 (Wyngaard) formula over water

			//jref:start commented out to see if stability is affected.
			double FLUXC = Math.Max(HFX / RHOX / Constants.cp + Constants.EP_1 * TSKV * QFX / RHOX, 0);
			double VCONVC = 1.0;
			double VCONV = VCONVC * Math.Pow(Constants.g / TGDSA * PBLH * FLUXC, 0.33);
			//  VCONV = 0
			//jref:end

			// Mahrt and Sun low-res correction
			double VSGD = 0.32 * Math.Pow(Math.Max(DX / 5000 - 1, 0), 0.33);
			WSPD = Math.Sqrt(WSPD * WSPD + VCONV * VCONV + VSGD * VSGD);
			if(double.IsNaN(WSPD))
				WSPD=0.1;
			WSPD = Math.Max(WSPD, 0.1);
			double BR = GOVRTH * ZA * DTHVDZ / (WSPD * WSPD);
			//  if previously unstable, do not let into regimes 1 and 2
			if (MOL < 0)
				BR = Math.Min(BR, 0.0);
			RMOL = -GOVRTH * DTHVDZ * ZA * Constants.KARMAN;
			//-----------------------------------------------------------------------
			//     diagnose basic parameters for the appropriated stability class:
//
			//     the stability classes are determined by br (bulk richardson no.)
			//     and hol (height of pbl/monin-obukhov length).
//
			//     criteria for the classes are as follows:
//
			//        1. br .ge. 0.2;
			//               represents nighttime stable conditions (regime=1),
//
			//        2. br < 0.2 && br .gt. 0.0;
			//               represents damped mechanical turbulent conditions
			//               (regime=2),
//
			//        3. br == 0.0
			//               represents forced convection conditions (regime=3),
//
			//        4. br < 0.0
			//               represents free convection conditions (regime=4).
//
			//!-----------------------------------------------------------------------
			int REGIME = 0;
			if (BR >= 0.2)
				REGIME = 1;
			if (BR < 0.2 && BR > 0.0)
				REGIME = 2;
			if (Math.Abs(BR)<1e-10)
				REGIME = 3;
			if (BR < 0.0)
				REGIME = 4;

			switch (REGIME) {
					
				case 1:
					{
						// class 1; stable (nighttime) conditions:
						PSIM = -10 * GZ1OZ0;
						// lower limit on psi in stable conditions
						PSIM = Math.Max(PSIM, -10.0);
						PSIH = PSIM;
						PSIM10 = 10.0 / ZA * PSIM;
						PSIM10 = Math.Max(PSIM10, -10.0);
						PSIH10 = PSIM10;
						double PSIM2 = 2 / ZA * PSIM;
						PSIM2 = Math.Max(PSIM2, -10.0);
						PSIH2 = PSIM2;

						// 1.0 over Monin-Obukhov length
						if (UST < 0.01)
							RMOL = BR * GZ1OZ0; //ZA/L
						else
							RMOL = Constants.KARMAN * GOVRTH * ZA * MOL / (UST * UST); //ZA/L
						
						RMOL = Math.Min(RMOL, 9.999); // ZA/L
						RMOL = RMOL / ZA; //1.0/L
						break;
					}
				case 2:
					{
						// class 2; damped mechanical turbulence:
						PSIM = -5.0 * BR * GZ1OZ0 / (1.1 - 5.0 * BR);
						// lower limit on psi in stable conditions
						PSIM = Math.Max(PSIM, -10.0);
						// AKB(1976), EQ(16).
						PSIH = PSIM;
						PSIM10 = 10.0 / ZA * PSIM;
						PSIM10 = Math.Max(PSIM10, -10.0);
						PSIH10 = PSIM10;
						double PSIM2 = 2 / ZA * PSIM;
						PSIM2 = Math.Max(PSIM2, -10.0);
						PSIH2 = PSIM2;

						// Linear form: PSIM = -0.5*ZA/L; e.g, see eqn 16 of
						// Blackadar, Modeling the nocturnal boundary layer, Preprints,
						// Third Symposium on Atmospheric Turbulence Diffusion and Air Quality,
						// Raleigh, NC, 1976
						double ZOL = BR * GZ1OZ0 / (1.00001 - 5.0 * BR);

						if (ZOL > 0.5) {  // linear form ok
							// Holtslag and de Bruin, J. App. Meteor 27, 689-704, 1988;
							// see also, Launiainen, Boundary-Layer Meteor 76,165-179, 1995
							// Eqn (8) of Launiainen, 1995
							ZOL = (1.89 * GZ1OZ0 + 44.2) * BR * BR + (1.18 * GZ1OZ0 - 1.37) * BR;
							ZOL = Math.Min(ZOL, 9.999);
						}
						// 1.0 over Monin-Obukhov length
						RMOL = ZOL / ZA;
						break;
					}

				case 3:  // class 3; forced convection:
					{
						PSIM = 0.0;
						PSIH = PSIM;
						PSIM10 = 0;
						PSIH10 = PSIM10;
						double	PSIM2 = 0;
						PSIH2 = PSIM2;
						double ZOL = Constants.KARMAN * GOVRTH * ZA * MOL / (UST * UST);
						if (UST < 0.01)
							ZOL = BR * GZ1OZ0;
						RMOL = ZOL / ZA;
						break;
					}
				case 4: // class 4; free convection:
					{
						double	ZOL = Constants.KARMAN * GOVRTH * ZA * MOL / (UST * UST);
						if (UST < 0.01)
							ZOL = BR * GZ1OZ0;
						
						
						double ZOL10 = 10.0 / ZA * ZOL;
						double ZOL2 = 2 / ZA * ZOL;
						ZOL = Math.Min(ZOL, 0);
						ZOL = Math.Max(ZOL, -9.9999);
						ZOL10 = Math.Min(ZOL10, 0);
						ZOL10 = Math.Max(ZOL10, -9.9999);
						ZOL2 = Math.Min(ZOL2, 0);
						ZOL2 = Math.Max(ZOL2, -9.9999);
						int NZOL = (int)(-ZOL * 100);
						double RZOL = -ZOL * 100 - NZOL;
						int NZOL10 = (int)(-ZOL10 * 100);
						double RZOL10 = -ZOL10 * 100 - NZOL10;
						int NZOL2 = (int)(-ZOL2 * 100);
						double RZOL2 = -ZOL2 * 100 - NZOL2;
						PSIM = Constants.PSIMTB[NZOL] + RZOL * (Constants.PSIMTB[NZOL + 1] - Constants.PSIMTB[NZOL]);
						PSIH = Constants.PSIHTB[NZOL] + RZOL * (Constants.PSIHTB[NZOL + 1] - Constants.PSIHTB[NZOL]);
						PSIM10 = Constants.PSIMTB[NZOL10] + RZOL10 * (Constants.PSIMTB[NZOL10 + 1] - Constants.PSIMTB[NZOL10]);
						PSIH10 = Constants.PSIHTB[NZOL10] + RZOL10 * (Constants.PSIHTB[NZOL10 + 1] - Constants.PSIHTB[NZOL10]);
						double PSIM2 = Constants.PSIMTB[NZOL2] + RZOL2 * (Constants.PSIMTB[NZOL2 + 1] - Constants.PSIMTB[NZOL2]);
						PSIH2 = Constants.PSIHTB[NZOL2] + RZOL2 * (Constants.PSIHTB[NZOL2 + 1] - Constants.PSIHTB[NZOL2]);

						// limit psih and psim in the case of thin layers and high roughness
						// this prevents denominator in fluxes from getting too small
						//       PSIH=Math.Min(PSIH,0.9*GZ1OZ0)
						//       PSIM=Math.Min(PSIM,0.9*GZ1OZ0)
						PSIH = Math.Min(PSIH, 0.9 * GZ1OZ0);
						PSIM = Math.Min(PSIM, 0.9 * GZ1OZ0);
						PSIH2 = Math.Min(PSIH2, 0.9 * GZ2OZ0);
						PSIM10 = Math.Min(PSIM10, 0.9 * GZ10OZ0);
						// AHW: mods to compute ck, cd
						PSIH10 = Math.Min(PSIH10, 0.9 * GZ10OZ0);

						RMOL = ZOL / ZA;
						break;
					}
			} // stability regime done
			// compute the frictional velocity: ZA(1982) EQS(2.60),(2.61).
			double DTG = THX - THGB;
			double PSIX = GZ1OZ0 - PSIM;
			double PSIX10 = GZ10OZ0 - PSIM10;

			// lower limit added to prevent large flhc in soil model
			// activates in unstable conditions with thin layers or high z0
			PSIT = Math.Max(GZ1OZ0 - PSIH, 2); //does this still apply???? jref
			double PSIQ = Math.Log(Constants.KARMAN * UST * ZA / XKA + ZA / ZL) - PSIH;
			double PSIT2 = GZ2OZ0 - PSIH2;
			double PSIQ2 = Math.Log(Constants.KARMAN * UST * 2 / XKA + 2 / ZL) - PSIH2;
			// AHW: mods to compute ck, cd
			double PSIQ10 = Math.Log(Constants.KARMAN * UST * 10.0 / XKA + 10.0 / ZL) - PSIH10;

			//jref:start - commented out since these values can be produced by sfclay routine
			//   if(PRESENT(ck) && PRESENT(cd) && PRESENT(cka) && PRESENT(cda))
			//      Ck=(karman/psix10)*(karman/psiq10)
			//      Cd=(karman/psix10)*(karman/psix10)
			//      Cka=(karman/psix)*(karman/psiq)
			//      Cda=(karman/psix)*(karman/psix)
			//   ENDIF

			//   WRITE(*,*) "KARMAN=",KARMAN
			//   WRITE(*,*) "UST=",UST
			//   WRITE(*,*) "XKA=",XKA
			//   WRITE(*,*) "ZA =",ZA
			//   WRITE(*,*) "ZL =",ZL
			//   WRITE(*,*) "PSIH=",PSIH
			//   WRITE(*,*) "PSIQ=",PSIQ,"PSIT=",PSIT
			
			if (IZ0TLND == 1) {
				ZL = ZNT;
				//       czil related changes for land
				double VISC = (1.32 + 0.009 * (T1D - 273.15)) * 1E-5;
				double RESTAR = UST * ZL / VISC;
				//       modify CZIL accordi;ng to Chen & Zhang, 2009

				double CZIL = Math.Pow(10.0, -0.40 * (ZL / 0.07));

				PSIT = GZ1OZ0 - PSIH + CZIL * Constants.KARMAN * Math.Sqrt(RESTAR);
				PSIQ = GZ1OZ0 - PSIH + CZIL * Constants.KARMAN * Math.Sqrt(RESTAR);
				PSIT2 = GZ2OZ0 - PSIH2 + CZIL * Constants.KARMAN * Math.Sqrt(RESTAR);
				PSIQ2 = GZ2OZ0 - PSIH2 + CZIL * Constants.KARMAN * Math.Sqrt(RESTAR);
			}


			// to prevent oscillations average with old value
			UST = 0.5 * UST + 0.5 * Constants.KARMAN * WSPD / PSIX;
//			if (double.IsNaN(UST))
//				throw new Exception();
			if(double.IsNaN(UST))
				UST=-1e20;
			UST = Math.Max(UST, 0.1);
			//jref: should this be converted to RMOL???
			MOL = Constants.KARMAN * DTG / PSIT / PRT;
			double DENOMQ = PSIQ;
			double DENOMQ2 = PSIQ2;
			double DENOMT2 = PSIT2;
			//   WRITE(*,*) "ILOC,JLOC=",ILOC,JLOC,"DENOMQ=",DENOMQ
			//   WRITE(*,*) "UST=",UST,"PSIT=",PSIT
			//   wrf_error_fatal("stop in sfcdif4")
			// calculate exchange coefficients
			//jref: start exchange coefficient for momentum
			CM = Constants.KARMAN * Constants.KARMAN / (PSIX * PSIX);
			//jref:end
			CHS = UST * Constants.KARMAN / DENOMQ;
			//        GZ2OZ0=ALOG(2./ZNT)
			//        PSIM2=-10.0*GZ2OZ0
			//        PSIM2=Math.Max(PSIM2,-10)
			//        PSIH2=PSIM2
			CQS2 = UST * Constants.KARMAN / DENOMQ2;
			CHS2 = UST * Constants.KARMAN / DENOMT2;
			// jref: in last iteration calculate diagnostics

			U10 = UX * PSIX10 / PSIX;
			V10 = VX * PSIX10 / PSIX;

			// jref: check the following for correct calculation
			//   TH2=THGB+DTG*PSIT2/PSIT
			//   Q2=QSFC+(QX-QSFC)*PSIQ2/PSIQ
			//   T2 = TH2*(PSFCPA/P1000mb)**RCP
			

		}
		
	}
}
