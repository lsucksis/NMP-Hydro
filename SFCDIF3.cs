/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2018/1/20
 * Time: 14:45
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;

namespace NoahMP
{
	/// <summary>
	/// Description of SFCDIF3.
	/// </summary>
	public static class SFCDIF3
	{
		public static FortDoubleArray PSIM2;
		public static FortDoubleArray PSIH2;
		
		static	double ZTMIN1 = -5.0;
		static		double ZTMAX1 = 1.0;
		static	double ZTMIN2 = -5.0;
		static		double ZTMAX2 = 1.0;
		static	double ZRNG1 = ZTMAX1 - ZTMIN1;
		static	double	ZRNG2 = ZTMAX2 - ZTMIN2;
		static int KZTM = 10001;
		static int KZTM2 = KZTM - 2;
		static	double DZETA1 = ZRNG1 / (KZTM - 1);
		static	double DZETA2 = ZRNG2 / (KZTM - 1);
		static	double ZETA1 = ZTMIN1;
		static	double ZETA2 = ZTMIN2;
		static	double P1000mb = 100000;
		
		public static void CalPSIM2(double ZETA2)
		{
			PSIH2 = new FortDoubleArray(1, KZTM);
			PSIM2 = new FortDoubleArray(1, KZTM);
			for (int K = 1; K <= KZTM; K++) {
				if (ZETA2 < 0) {
//
//----------------------------------------------------------------------
//***  PAULSON 1970 FUNCTIONS
//----------------------------------------------------------------------
//
					double X = Math.Sqrt(Math.Sqrt(1 - 16 * ZETA2));
//
					PSIM2[K] = -2 * Math.Log((X + 1) / 2) - Math.Log((X * X + 1) / 2) + 2 * Math.Atan(X) - Math.PI / 2;
					PSIH2[K] = -2 * Math.Log((X * X + 1) / 2);
//----------------------------------------------------------------------
//***  STABLE RANGE
//----------------------------------------------------------------------
//
				} else {
//
//----------------------------------------------------------------------
//***  PAULSON 1970 FUNCTIONS
//----------------------------------------------------------------------
//
//         PSIM2(K)=5.*ZETA2
//         PSIH2(K)=5.*ZETA2
//
//----------------------------------------------------------------------
//***  HOLTSLAG AND DE BRUIN 1988
//----------------------------------------------------------------------
//
					PSIM2[K] = 0.7 * ZETA2 + 0.75 * ZETA2 * (6 - 0.35 * ZETA2) * Math.Exp(-0.35 * ZETA2);
					PSIH2[K] = 0.7 * ZETA2 + 0.75 * ZETA2 * (6 - 0.35 * ZETA2) * Math.Exp(-0.35 * ZETA2);
//----------------------------------------------------------------------
//
				}

			}
		}
		/// <summary>
		///   computing surface drag coefficient CM for momentum and CH for heat
		///  Joakim Refslund, 2011, MYJ SFCLAY
		/// </summary>
		/// <param name="ILOC"></param>
		/// <param name="JLOC"></param>
		/// <param name="TSK"></param>
		/// <param name="QS"></param>
		/// <param name=""></param>
		/// <param name="PBLH"></param>
		/// <param name="Z0"></param>
		/// <param name="Z0BASE"></param>
		/// <param name="VEGTYP"></param>
		/// <param name="ISURBAN"></param>
		/// <param name="IZ0TLND"></param>
		/// <param name="SFCSPD"></param>
		/// <param name="ITER"></param>
		/// <param name="ITRMX"></param>
		/// <param name="TLOW"></param>
		/// <param name="THLOW"></param>
		/// <param name="QLOW"></param>
		/// <param name="CWMLOW"></param>
		/// <param name="ZSL"></param>
		/// <param name="PLOW"></param>
		/// <param name="USTAR"></param>
		/// <param name="AKMS"></param>
		/// <param name="AKHS"></param>
		/// <param name="CHS2"></param>
		/// <param name="CQS2"></param>
		/// <param name="RLMO"></param>
		/// <returns></returns>
		public static void SFCDIF(int ILOC, int JLOC, double TSK, double QS, double PSFC,  //in  
			double   PBLH, double Z0, double Z0BASE, double VEGTYP, int ISURBAN,  //in  
			int  IZ0TLND, double SFCSPD, int ITER, int ITRMX, double TLOW,  //in  
			double  THLOW, double  QLOW, double CWMLOW, double ZSL,          //in
			double   PLOW, ref double USTAR, ref double AKMS, ref double AKHS, out double CHS2,  //inout
			out double  CQS2, out double RLMO)                           //out 
		{                     

			double EPSU2 = 1E-6;
			double EPSUST = 1e-9; 
			double EPSZT = 1E-28;
			double A2S = 17.2693882, A3S = 273.16, A4S = 35.86;
			double SEAFC = 0.98;
			//double PQ0SEA=PQ0*SEAFC;
			double BETA = 1 / 273, EXCML = 0.0001, EXCMS = 0.0001;
			double GLKBR = 10, GLKBS = 30, PI = 3.1415926;
			double QVISC = 2.1E-5, RIC = 0.505, SMALL = 0.35;
			double SQPR = 0.84, SQSC = 0.84, SQVISC = 258.2, TVISC = 2.1E-5;
			double USTC = 0.7, USTR = 0.225, VISC = 1.5E-5, FH = 1.01;
			double WWST = 1.2, ZTFC = 1, TOPOFAC = 9.0e-6;
                

//----------------------------------------------------------------------   

//----------------------------------------------------------------------   

//    INTEGER,INTENT(IN) :: ILOC
//    INTEGER,INTENT(IN) :: JLOC
//    REAL   ,INTENT(IN) :: TSK
//    REAL   ,INTENT(IN) :: PSFC
//    REAL   ,INTENT(IN) :: PBLH
//    INTEGER,INTENT(IN) :: VEGTYP  //in routine
//    INTEGER,INTENT(IN) :: ISURBAN //in veg_parm                    
//    INTEGER,INTENT(IN) :: IZ0TLND
//    REAL   ,INTENT(IN) :: QLOW
//    REAL   ,INTENT(IN) :: THLOW
//    REAL   ,INTENT(IN) :: TLOW
//    REAL   ,INTENT(IN) :: CWMLOW
//    REAL   ,INTENT(IN) :: SFCSPD
//    REAL   ,INTENT(IN) :: PLOW
//    REAL   ,INTENT(IN) :: ZSL
//    REAL   ,INTENT(IN) :: Z0BASE
//    INTEGER,INTENT(IN) :: ITER
//    INTEGER,INTENT(IN) :: ITRMX
//
//// output                                                                        
//    REAL   ,INTENT(OUT) :: CHS2
//    REAL   ,INTENT(OUT) :: CQS2
//    REAL   ,INTENT(OUT) :: RLMO
//
//// input/output
//    REAL   ,INTENT(INOUT) :: AKHS
//    REAL   ,INTENT(INOUT) :: AKMS
//    REAL   :: QZ0
//    REAL   ,INTENT(INOUT) :: USTAR
//    REAL   ,INTENT(IN) :: Z0
//    REAL   ,INTENT(INOUT):: QS
//    REAL                :: RIB

// local
//    INTEGER :: ITR,K
//    REAL :: THZ0
//    REAL :: THVLOW
//    REAL :: CT
//    REAL :: BTGH
//    REAL :: BTGX
//    REAL :: CXCHL
//    REAL :: DTHV
//    REAL :: DU2
//    REAL :: ELFC
//    REAL :: PSH02
//    REAL :: PSH10
//    REAL :: PSHZ
//    REAL :: PSHZL
//    REAL :: PSM10
//    REAL :: PSMZ
//    REAL :: PSMZL
//    REAL :: RDZ
//    REAL :: RDZT
//			   REAL :: RLMA //???
//    REAL :: RLMN //???
//    REAL :: RLOGT
//    REAL :: RLOGU
//    REAL :: RZ
//    REAL :: SIMH
//    REAL :: SIMM
//    REAL :: USTARK
//    REAL :: WSTAR2
//    REAL :: WSTAR
//    REAL :: CHS
//    REAL :: RZSU
//    REAL :: RZST
//    REAL :: X,XLT,XLT4,XLU,XLU4,XT,XT4,XU,XU4,ZETALT,ZETALU         , 
//            ZETAT,ZETAU,ZQ,ZSLT,ZSLU,ZT,ZU,TOPOTERM,ZZIL
//    REAL :: AKHS02,AKHS10,AKMS02,AKMS10
//    REAL :: ZU10
//    REAL :: ZT02
//    REAL :: ZT10
//    REAL :: RLNU10
//    REAL :: RLNT02
//    REAL :: RLNT10
//    REAL :: ZTAU10
//    REAL :: ZTAT02
//    REAL :: ZTAT10
//    REAL :: SIMM10
//    REAL :: SIMH02
//    REAL :: SIMH10
//    REAL :: ZUUZ
//    REAL :: EKMS10
//    REAL :: test
//    REAL :: E1
			CalPSIM2(ZETA2);
			double VKRM = 0.40;
			double CZETMAX = 10;

// diagnostic terms                                                         

//   REAL :: CZIL
//   REAL :: ZILFC

//double KTMZ,KTMZ2,DZETA1,DZETA2,FH01,FH02,ZTMAX1,ZTMAX2
	
// PSIH1,PSIH2,PSIM1,PSIM2 ARE DEFINED IN MODULE_SF_MYJSFC
			

// calculate potential and virtual potential temperatures
			double THVLOW = THLOW * (1 + Constants.EP_1 * QLOW);
			double THZ0 = TSK * Math.Pow(P1000mb / PSFC, Constants.rcp);

// calculate initial values
			double ZU = Z0;
			double ZT = ZU * ZTFC;      //ZTFC = ZOH/ZOM =<1 set to 1 at beginning
			double ZQ = ZT;
			double QZ0 = QS;

			double RDZ = 1 / ZSL;
			double CXCHL = EXCML * RDZ;
			double DTHV = THVLOW - THZ0 * (0.608 * QZ0 + 1);    //delta pot. virtual temperature

			double BTGX = NoahMP.GRAV / THLOW;
			double ELFC = VKRM * BTGX;
			double BTGH = 0;
// Minimum PBLH is >= 1000.
			if (PBLH > 1000) {
				BTGH = BTGX * PBLH;
			} else {
				BTGH = BTGX * 1000;
			}

			double DU2 = Math.Max(SFCSPD * SFCSPD, EPSU2);  //Wind speed - EPSU2 parm = 1*10^-6
			double RIB = BTGX * DTHV * ZSL / DU2;         //Bulk richardson stability
			double ZSLU = ZSL + ZU;
			double RZSU = ZSLU / ZU;
			double RLOGU = Math.Log(RZSU);       //log(z/z0)

			double ZSLT = ZSL + ZU;
			double CZIL = 0;
			if ((IZ0TLND == 0) || (VEGTYP == ISURBAN)) {    // ARE IZ0TLND DEFINED HERE?           
				// Just use the original CZIL value.                          
				CZIL = 0.1;
			} else {
				// Modify CZIL according to Chen  Zhang, 2009                
				// CZIL = 10 ** -0.40 H, ( where H = 10*Zo )                  
				CZIL = Math.Pow(10.0, -0.40 * (Z0 / 0.07));
			}
			double	ZILFC = -CZIL * VKRM * SQVISC;     //SQVISC parm

			double ZZIL = 0;
			// stable   				                                                               
			if (DTHV > 0) {
				if (RIB < RIC)
					ZZIL = ZILFC * (1.0 + (RIB / RIC) * (RIB / RIC) * CZETMAX);
				else
					ZZIL = ZILFC * (1.0 + CZETMAX);
			} else {// unstable
				ZZIL = ZILFC;
			}

//---  ZILITINKEVITCH FIX FOR ZT                                           
// oldform   ZT=Math.Max(EXP(ZZIL*Math.Sqrt(USTAR*ZU))*ZU,EPSZT)                     
			ZT = Math.Max(Math.Exp(ZZIL * Math.Sqrt(USTAR * Z0BASE)) * Z0BASE, EPSZT);  //Z0 is backgrund roughness?
			double RZST = ZSLT / ZT;
			double RLOGT = Math.Log(RZST);
	
//----------------------------------------------------------------------   
//  1./MONIN-OBUKHOV LENGTH-SCALE                                       
//----------------------------------------------------------------------   
			RLMO = ELFC * AKHS * DTHV / (USTAR * USTAR * USTAR);

			double ZETALU = ZSLU * RLMO;
			double ZETALT = ZSLT * RLMO;
			double ZETAU = ZU * RLMO;
			double ZETAT = ZT * RLMO;

			ZETALU = Math.Min(Math.Max(ZETALU, ZTMIN2), ZTMAX2);
			ZETALT = Math.Min(Math.Max(ZETALT, ZTMIN2), ZTMAX2);
			ZETAU = Math.Min(Math.Max(ZETAU, ZTMIN2 / RZSU), ZTMAX2 / RZSU);
			ZETAT = Math.Min(Math.Max(ZETAT, ZTMIN2 / RZST), ZTMAX2 / RZST);

//----------------------------------------------------------------------   
//***  LAND FUNCTIONS                                                      
//----------------------------------------------------------------------   

			double RZ = (ZETAU - ZTMIN2) / DZETA2;
			int K = (int)(RZ);
			double RDZT = RZ - (double)(K);
			K = Math.Min(K, KZTM2);
			K = Math.Max(K, 0);
			double	PSMZ = (PSIM2[K + 2] - PSIM2[K + 1]) * RDZT + PSIM2[K + 1];

			RZ = (ZETALU - ZTMIN2) / DZETA2;
			K = (int)(RZ);
			RDZT = RZ - (double)(K);
			K = Math.Min(K, KZTM2);
			K = Math.Max(K, 0);
			double PSMZL = (PSIM2[K + 2] - PSIM2[K + 1]) * RDZT + PSIM2[K + 1];

			double SIMM = PSMZL - PSMZ + RLOGU;

			RZ = (ZETAT - ZTMIN2) / DZETA2;
			K = (int)(RZ);
			RDZT = RZ - (double)(K);
			K = Math.Min(K, KZTM2);
			K = Math.Max(K, 0);
			double	PSHZ = (PSIH2[K + 2] - PSIH2[K + 1]) * RDZT + PSIH2[K + 1];

			RZ = (ZETALT - ZTMIN2) / DZETA2;
			K = (int)(RZ);
			RDZT = RZ - (double)(K);
			K = Math.Min(K, KZTM2);
			K = Math.Max(K, 0);
			double PSHZL = (PSIH2[K + 2] - PSIH2[K + 1]) * RDZT + PSIH2[K + 1];
			double FH02 = 1.0;
			double SIMH = (PSHZL - PSHZ + RLOGT) * FH02;
//----------------------------------------------------------------------   
			double USTARK = USTAR * VKRM;
			AKMS = Math.Max(USTARK / SIMM, CXCHL);
			AKHS = Math.Max(USTARK / SIMH, CXCHL);
			//----------------------------------------------------------------------   
//  BELJAARS CORRECTION FOR USTAR                                       
//----------------------------------------------------------------------   
			double WSTAR2 = 0;
			if (DTHV <= 0) {
				double WWST2 = WWST * WWST;
				WSTAR2 = WWST2 * Math.Pow(Math.Abs(BTGH * AKHS * DTHV), 2.0 / 3);                //zj
			} else {
				WSTAR2 = 0;                                              
			}
			//zj

			USTAR = Math.Max(Math.Sqrt(AKMS * Math.Sqrt(DU2 + WSTAR2)), EPSUST);
			double CT = 0;
//----------------------------------------------------------------------   
//***  THE FOLLOWING DIAGNOSTIC BLOCK PRODUCES 2-m and 10-m VALUES         
//***  FOR TEMPERATURE, MOISTURE, AND WINDS.  IT IS DONE HERE SINCE        
//***  THE VARIOUS QUANTITIES NEEDED FOR THE COMPUTATION ARE LOST          
//***  UPON EXIT FROM THE ROTUINE.                                         
//----------------------------------------------------------------------   
			double WSTAR = Math.Sqrt(WSTAR2) / WWST;

//jref: calculate in last iteration
//  if (ITER == ITRMX) {

			double ZU10 = ZU + 10;
			double ZT02 = ZT + 02;
			double ZT10 = ZT + 10;

			double RLNU10 = Math.Log(ZU10 / ZU);
			double RLNT02 = Math.Log(ZT02 / ZT);
			double RLNT10 = Math.Log(ZT10 / ZT);

			double ZTAU10 = ZU10 * RLMO;
			double ZTAT02 = ZT02 * RLMO;
			double ZTAT10 = ZT10 * RLMO;

			ZTAU10 = Math.Min(Math.Max(ZTAU10, ZTMIN2), ZTMAX2);
			ZTAT02 = Math.Min(Math.Max(ZTAT02, ZTMIN2), ZTMAX2);
			ZTAT10 = Math.Min(Math.Max(ZTAT10, ZTMIN2), ZTMAX2);

//jref: land diagnostic functions
			RZ = (ZTAU10 - ZTMIN2) / DZETA2;
			K = (int)(RZ);
			RDZT = RZ - (double)(K);
			K = Math.Min(K, KZTM2);
			K = Math.Max(K, 0);
			double PSM10 = (PSIM2[K + 2] - PSIM2[K + 1]) * RDZT + PSIM2[K + 1];

			double SIMM10 = PSM10 - PSMZ + RLNU10;

			RZ = (ZTAT02 - ZTMIN2) / DZETA2;
			K = (int)(RZ);
			RDZT = RZ - (double)(K);
			K = Math.Min(K, KZTM2);
			K = Math.Max(K, 0);
			double PSH02 = (PSIH2[K + 2] - PSIH2[K + 1]) * RDZT + PSIH2[K + 1];

			double SIMH02 = (PSH02 - PSHZ + RLNT02) * FH02;

			RZ = (ZTAT10 - ZTMIN2) / DZETA2;
			K = (int)(RZ);
			RDZT = RZ - (double)(K);
			K = Math.Min(K, KZTM2);
			K = Math.Max(K, 0);
			double PSH10 = (PSIH2[K + 2] - PSIH2[K + 1]) * RDZT + PSIH2[K + 1];

			double SIMH10 = (PSH10 - PSHZ + RLNT10) * FH02;
//jref: diagnostic exchange coefficients                                                                    
			double AKMS10 = Math.Max(USTARK / SIMM10, CXCHL);
			double AKHS02 = Math.Max(USTARK / SIMH02, CXCHL);
			double AKHS10 = Math.Max(USTARK / SIMH10, CXCHL);

			double ZUUZ = Math.Min(ZU * 0.50, 0.18);
			ZU = Math.Max(ZU * 0.35, ZUUZ);
//                                                                         
			ZU10 = ZU + 10;
			RZSU = ZU10 / ZU;
			RLNU10 = Math.Log(RZSU);

			ZETAU = ZU * RLMO;
			ZTAU10 = ZU10 * RLMO;

			ZTAU10 = Math.Min(Math.Max(ZTAU10, ZTMIN2), ZTMAX2);
			ZETAU = Math.Min(Math.Max(ZETAU, ZTMIN2 / RZSU), ZTMAX2 / RZSU);

			RZ = (ZTAU10 - ZTMIN2) / DZETA2;
			K = (int)(RZ);
			RDZT = RZ - (double)(K);
			K = Math.Min(K, KZTM2);
			K = Math.Max(K, 0);
			PSM10 = (PSIM2[K + 2] - PSIM2[K + 1]) * RDZT + PSIM2[K + 1];
			SIMM10 = PSM10 - PSMZ + RLNU10;
			double	EKMS10 = Math.Max(USTARK / SIMM10, CXCHL);

//        U10E=UMFLX/EKMS10+UZ0                                             
//        V10E=VMFLX/EKMS10+VZ0                                             

//      ENDif                                                               
//                                                                         
//      U10=U10E                                                            
//      V10=V10E                                                            
//                                                                         
//----------------------------------------------------------------------   
//***  SET OTHER WRF DRIVER ARRAYS                                         
//----------------------------------------------------------------------  
			//                                                                         
//jref commented out
//      RLOW=PLOW/(R_D*TLOW)                                                
			double	CHS = AKHS;
			CHS2 = AKHS02;
			CQS2 = AKHS02;

//  END if
		}
	}
}