/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2016/11/27
 * Time: 8:11
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;

namespace NoahMP
{
	/// <summary>
	/// Description of VEG_PARAMS.
	/// </summary>
	public class VEG_PARAMS
	{
		int MAX_VEG_PARAMS = 33;
		/// <summary>
		/// No. of vegetation types
		/// </summary>
		public static int MVT = 27;
		public static int MBAND = 2;

		public int ISURBAN;
		public int ISWATER;
		public int ISBARREN;
		public int ISSNOW;
		public int EBLFOREST;

		/// <summary>
		/// maximum intercepted water per unit lai+sai (mm)
		/// </summary>
		public double[] CH2OP = new double[MVT];
		/// <summary>
		/// characteristic leaf dimension (m)
		/// </summary>
		public double[] DLEAF = new double[MVT];
		/// <summary>
		/// momentum roughness length (m)
		/// </summary>
		public double[] Z0MVT = new double[MVT];
		/// <summary>
		/// top of canopy (m)
		/// </summary>
		public double[] HVT = new double[MVT];
		/// <summary>
		/// bottom of canopy (m)
		/// </summary>
		public double[] HVB = new double[MVT];
		/// <summary>
		/// tree density (no. of trunks per m2)
		/// </summary>
		public double[] DEN = new double[MVT];
		/// <summary>
		/// tree crown radius (m)
		/// </summary>
		public double[] RC = new double[MVT];
		/// <summary>
		/// monthly stem area index, one-sided
		/// </summary>
		public double[,] SAIM = new double[MVT, 12];
		/// <summary>
		/// monthly leaf area index, one-sided
		/// </summary>
		public double[,] LAIM = new double[MVT, 12];
		/// <summary>
		/// single-side leaf area per Kg [m2/kg]
		/// </summary>
		public double[] SLA = new double[MVT];
		/// <summary>
		/// coeficient for leaf stress death [1/s]
		/// </summary>
		public double[] DILEFC = new double[MVT];
		/// <summary>
		/// coeficient for leaf stress death [1/s]
		/// </summary>
		public double[] DILEFW = new double[MVT];
		/// <summary>
		/// fraction of growth respiration  //original was 0.3
		/// </summary>
		public double[] FRAGR = new double[MVT];
		/// <summary>
		/// leaf turnover [1/s]
		/// </summary>
		public double[] LTOVRC = new double[MVT];
		/// <summary>
		/// photosynthetic pathway: 0. = c4, 1. = c3
		/// </summary>
		public double[] C3PSN = new double[MVT];
		/// <summary>
		/// co2 michaelis-menten constant at 25c (pa)
		/// </summary>
		public double[] KC25 = new double[MVT];
		/// <summary>
		/// q10 for kc25
		/// </summary>
		public double[] AKC = new double[MVT];
		/// <summary>
		/// o2 michaelis-menten constant at 25c (pa)
		/// </summary>
		public double[] KO25 = new double[MVT];
		/// <summary>
		/// q10 for ko25
		/// </summary>
		public double[] AKO = new double[MVT];
		
		/// <summary>
		/// maximum rate of carboxylation at 25c (umol co2/m**2/s)
		/// </summary>
		public double[] VCMX25 = new double[MVT];
		
		/// <summary>
		/// q10 for vcmx25
		/// </summary>
		public double[] AVCMX = new double[MVT];
		
		/// <summary>
		/// minimum leaf conductance (umol/m**2/s)
		/// </summary>
		public double[] BP = new double[MVT];
		
		/// <summary>
		/// slope of conductance-to-photosynthesis relationship
		/// </summary>
		public double[] MP = new double[MVT];
		
		/// <summary>
		/// quantum efficiency at 25c (umol co2 / umol photon)
		/// </summary>
		public double[] QE25 = new double[MVT];
		
		/// <summary>
		/// q10 for qe25
		/// </summary>
		public double[] AQE = new double[MVT];
		
		/// <summary>
		/// leaf maintenance respiration at 25c (umol co2/m**2/s)
		/// </summary>
		public double[] RMF25 = new double[MVT];
		
		/// <summary>
		/// stem maintenance respiration at 25c (umol co2/kg bio/s)
		/// </summary>
		public double[] RMS25 = new double[MVT];
		
		/// <summary>
		/// root maintenance respiration at 25c (umol co2/kg bio/s)
		/// </summary>
		public double[] RMR25 = new double[MVT];
		
		/// <summary>
		/// q10 for maintenance respiration
		/// </summary>
		public double[] ARM = new double[MVT];
		
		/// <summary>
		/// foliage nitrogen concentration when f(n)=1 (%)
		/// </summary>
		public double[] FOLNMX = new double[MVT];
		
		/// <summary>
		/// minimum temperature for photosynthesis [K]
		/// </summary>
		public double[] TMIN = new double[MVT];
		/// <summary>
		/// leaf/stem orientation index
		/// </summary>

		public double[] XL = new double[MVT];
		/// <summary>
		/// leaf reflectance: 1=vis, 2=nir
		/// </summary>
		public double[,] RHOL = new double[MVT, MBAND];
		/// <summary>
		/// stem reflectance: 1=vis, 2=nir
		/// </summary>
		public double[,] RHOS = new double[MVT, MBAND];
		/// <summary>
		/// leaf transmittance: 1=vis, 2=nir
		/// </summary>
		public double[,] TAUL = new double[MVT, MBAND];
		/// <summary>
		/// stem transmittance: 1=vis, 2=nir
		/// </summary>
		public double[,] TAUS = new double[MVT, MBAND];
		/// <summary>
		/// microbial respiration parameter (umol co2 /kg c/ s)
		/// </summary>
		public double[] MRP = new double[MVT];
		
		/// <summary>
		/// empirical canopy wind parameter
		/// </summary>
		public double[] CWPVT = new double[MVT];
		/// <summary>
		/// wood to non-wood ratio
		/// </summary>
		public double[] WRRAT = new double[MVT];
		/// <summary>
		/// wood pool (switch 1 or 0) depending on woody or not [-]
		/// </summary>
		public double[] WDPOOL = new double[MVT];
		/// <summary>
		/// characteristic T for leaf freezing [K]
		/// </summary>
		public double[] TDLEF = new double[MVT];
		

		public int IK = 0;
		public int IM = 0;
		

		public double[] slarea = new double[MVT];
		public double[,] eps =null;// new double[MVT, 5];

		public VEG_PARAMS()
		{
		}
	}
}
