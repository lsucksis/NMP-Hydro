/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2018/1/12
 * Time: 14:15
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;

namespace NoahMP
{
	/// <summary>
	/// 这些参数需要从参数表.tbl中读取
	/// </summary>
	public class REDPRM
	{
		public REDPRM()
		{
		}
		/// <summary>
		/// maximum stomatal resistance ( Assigned in REDPRM )
		/// </summary>
		public static double RSMAX = 5000;
		/// <summary>
		/// minimum Canopy Resistance [s/m] ( Assigned in REDPRM )
		/// </summary>
		public static double RSMIN = 0;
		/// <summary>
		/// parameter used in vapor pressure deficit function ( Assigned in REDPRM )
		/// </summary>
		public static double HS = 0;
		/// <summary>
		/// parameter used in radiation stress function ( Assigned in REDPRM )
		/// </summary>
		public static double RGL = 0;
		/// <summary>
		/// optimum transpiration air temperature.
		/// </summary>
		public static double TOPT = 298;
		
		public static double[] RGL_Data;
		
		
		//		  INTEGER :: NROOT        //rooting depth [as the number of layers] ( Assigned in REDPRM )
		

		public static string LUTYPE;
		public static int LUCATS, BARE, NATURAL;
		public static int NLUS = 50;
		public static int[] NROTBL;
		public static double[] RSTBL, RGLTBL, HSTBL;
		public static double TOPT_DATA, RSMAX_DATA;

		//// not further used in this version (niu):

		public static double[] SNUPTBL, LAITBL, ALBTBL, SHDTBL, MAXALB;
		public static double CMCMAX_DATA, CFACTR_DATA, SBETA_DATA, SALP_DATA, SMLOW_DATA, SMHIGH_DATA, LVCOEF_DATA;

		public static	double[] LAIMINTBL;
		//KWM
		public static	double[]	LAIMAXTBL;
		//KWM
		public static	double[] EMISSMINTBL;
		//KWM
		public static double[] EMISSMAXTBL;
		//KWM
		public static double[] ALBEDOMINTBL;
		//KWM
		public static double[] ALBEDOMAXTBL;
		//KWM
		public static double[] Z0MINTBL;
		//KWM
		public static double[] Z0MAXTBL;
		//KWM
		public static double[] ZTOPVTBL;
		public static double[] ZBOTVTBL;
		
		public static double[] SLOPE_DATA;
//		public static double SBETA_DATA;
		public static double FXEXP_DATA;
		public static double CSOIL_DATA;
//		public static double SALP_DATA;
		public static double REFDK_DATA;
		public static double REFKDT_DATA;
		public static double FRZK_DATA;
		public static double ZBOT_DATA;
		public static double CZIL_DATA;
//		public static double SMLOW_DATA;
//		public static double SMHIGH_DATA;
//		public static double LVCOEF_DATA;
		public static void InitData(int n)
		{
			LUCATS = n;
			REDPRM.SHDTBL = new double[n];
			REDPRM.NROTBL = new int[n];
			REDPRM.RSTBL = new double[n];
			REDPRM.RGLTBL = new double[n];
			REDPRM.HSTBL = new double[n];
			REDPRM.SNUPTBL = new double[n];
			REDPRM.MAXALB = new double[n];
			REDPRM.LAIMINTBL = new double[n];
			REDPRM.LAIMAXTBL = new double[n];
			REDPRM.EMISSMINTBL = new double[n];
			REDPRM.EMISSMAXTBL = new double[n];
			REDPRM.ALBEDOMINTBL = new double[n];
			REDPRM.ALBEDOMAXTBL = new double[n];
			REDPRM.Z0MINTBL = new double[n];
			REDPRM.Z0MAXTBL = new double[n];
			REDPRM.ZTOPVTBL = new double[n];
			REDPRM.ZBOTVTBL = new double[n];
		}
	}
}
