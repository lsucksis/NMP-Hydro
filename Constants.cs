/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2018/1/12
 * Time: 9:23
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;

namespace NoahMP
{
	/// <summary>
	/// Description of Constants.
	/// </summary>
	public class Constants
	{
		
		//  2. Following are constants for use in defining real number bounds.

		//  A really small number.

		public static double epsilon = 1E-15;

		//  4. Following is information related to the physical constants.

		//  These are the physical constants used within the model.

		// JM NOTE -- can we name this grav instead?
		public static double g = 9.81;
		// acceleration due to gravity (m {s}^-2)


		public static double R_d = 287.04;
		public static double cp = 1004.6;

		//public static double  R_d          = 287;
		//public static double  cp           = 7*R_d/2;


		public static double R_v = 461.6;
		public static double cv = cp - R_d;
		public static double cpv = 4 * R_v;
		public static double cvv = cpv - R_v;
		public static double cvpm = -cv / cp;
		public static double cliq = 4190;
		public static double cice = 2106;
		public static double psat = 610.78;
		public static double rcv = R_d / cv;
		public static double rcp = R_d / cp;
		public static double rovg = R_d / g;
		public static double c2 = cp * rcv;
		public static double mwdry = 28.966;
		// molecular weight of dry air (g/mole)

		public static double p1000mb = 100000;
		public static double t0 = 300;
		public static double p0 = p1000mb;
		public static double cpovcv = cp / (cp - R_d);
		public static double cvovcp = 1 / cpovcv;
		public static double rvovrd = R_v / R_d;

		public static double reradius = 1 / 6370.0e03;

		public static double asselin = 0.025;
		//   public static double  asselin      = .0
		public static double cb = 25;

		public static double XLV0 = 3.15E6;
		public static double XLV1 = 2370;
		public static double XLS0 = 2.905E6;
		public static double XLS1 = 259.532;

		public static double XLS = 2.85E6;
		public static double XLV = 2.5E6;
		public static double XLF = 3.50E5;

		public static double rhowater = 1000;
		public static double rhosnow = 100;
		public static double rhoair0 = 1.28;
		//
		public static double n_ccn0 = 1.0E8;
		//
		public static double piconst = 3.1415926535897932384626433;
		public static double DEGRAD = piconst / 180;
		public static double DPD = 360 / 365;

		public static double SVP1 = 0.6112;
		public static double SVP2 = 17.67;
		public static double SVP3 = 29.65;
		public static double SVPT0 = 273.15;
		public static double   EP_1 { 
			get {
				return R_v / R_d - 1;
			}
		}
		public static double EP_2 = R_d / R_v;
		public static double KARMAN = 0.4;
		public static double EOMEG = 7.2921E-5;
		public static double STBOLT = 5.67051E-8;

		public static double prandtl = 1 / 3.0;
		// constants for w-damping option
		public static double w_alpha = 0.3;
		// strength m/s/s
		public static double w_beta = 1.0;
		// activation cfl number

		public static double pq0 = 379.90516;
		public static double epsq2 = 0.2;
		public static double a2 = 17.2693882;
		public static double a3 = 273.16;
		public static double a4 = 35.86;
		public static double epsq = 1E-12;
		public static double p608 = rvovrd - 1;
		//#if ( NMM_CORE == 1 )
		public static double climit = 1E-20;
		public static double cm1 = 2937.4;
		public static double cm2 = 4.9283;
		public static double cm3 = 23.5518;
		//       public static double   defc=8.0
		//       public static double   defm=32.0
		public static double defc = 0.0;
		public static double defm = 99999.0;
		public static double epsfc = 1 / 1.05;
		public static double epswet = 0.0;
		public static double fcdif = 1 / 3.0;
		//#ifdef HWRF
		public static double fcm = 0.0;
		//#else
		//public static double   fcm=0.00003;
		//#endif
		public static double gma = -R_d * (1 - rcp) * 0.5;
		public static double p400 = 40000.0;
		public static double phitp = 15000.0;
		public static double pi2 = 2 * 3.1415926;
		public static double	pi1 = 3.1415926;
		public static double plbtm = 105000.0;
		public static double plomd = 64200.0;
		public static double pmdhi = 35000.0;
		public static double q2ini = 0.50;
		public static double rfcp = 0.25 / cp;
		public static double rhcrit_land = 0.75;
		public static double rhcrit_sea = 0.80;
		public static double rlag = 14.8125;
		public static double rlx = 0.90;
		public static double scq2 = 50.0;
		public static double slopht = 0.001;
		public static double tlc = 2 * 0.703972477;
		public static double wa = 0.15;
		public static double wght = 0.35;
		public static double wpc = 0.075;
		public static double z0land = 0.10;
		//#ifdef HWRF
		public static double z0max = 0.01;
		//#else
		//public static double   z0max=0.008;
		//#endif
		public static double z0sea = 0.001;
		//#endif


		//  Earth

		//  The value for P2SI *must* be set to 1.0 for Earth
		//  Although, now we may not need this declaration here (see above)
		//public static double  P2SI         = 1.0

		//  Orbital constants:

		int PLANET_YEAR = 365;
		public static double OBLIQUITY = 23.5;
		public static double ECCENTRICITY = 0.014;
		public static double SEMIMAJORAXIS = 1.0;
		// In AU
		// Don't know the following values, so we'll fake them for now
		public static double zero_date = 0.0;
		// Time of perihelion passage
		//  Fraction into the year (from perhelion) of the
		//  occurrence of the Northern Spring Equinox
		public static double EQUINOX_FRACTION = 0.0;

		// 2012103
		//#if (EM_CORE == 1)
		// for calls to set_tiles
		int ZONE_SOLVE_EM = 1;
		int ZONE_SFS = 2;
		//#endif
		
		public static double[] PSIMTB=new double[1000];
		public static double[] PSIHTB=new double[1000];
		public Constants()
		{
			for(int i=0;i<1000;i++)
			{
				double ZOLN=-i*0.01;
				double X=Math.Pow(1-16*ZOLN,0.25);
				PSIMTB[i]=2*Math.Log(0.5*(1+X))+Math.Log(0.5*(1+X*X))-2*Math.Atan(X)+2*Math.Atan(1);
				double Y=Math.Pow(1-16*ZOLN,0.5);
				PSIHTB[i]=2*Math.Log(0.5*(1+Y));
			}
		}

	}
}
