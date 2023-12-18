/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2018/1/10
 * Time: 17:30
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;

namespace NoahMP
{
	/// <summary>
	/// Description of RAD_PARAMS.
	/// </summary>
	public class RAD_PARAMS
	{
		public static int MSC   = 9;
		public static int MBAND = 2;

		public static double[,] ALBSAT=new double[,]{{0.15,0.3},{0.11,0.22},{0.10,0.2},{0.09,0.18},{0.08,0.16},{0.07,0.14},{0.06,0.12},{0.05,0.1}};   //saturated soil albedos: 1=vis, 2=nir
		public static double[,] ALBDRY=new double[,]{{0.27,0.54},{0.22,0.44},{0.20,0.4},{0.18,0.36},{0.16,0.32},{0.14,0.28},{0.12,0.24},{0.10,0.2}};   //dry soil albedos: 1=vis, 2=nir
		public static double[] ALBICE=new double[]{0.8,0.55};       //albedo land ice: 1=vis, 2=nir
		public static double[] ALBLAK=new double[]{0.6,0.4};       //albedo frozen lakes: 1=vis, 2=nir
		public static double[] OMEGAS=new double[]{0.8,0.4};       //two-stream parameter omega for snow
		public static double BETADS=0.5;              //two-stream parameter betad for snow
		public static double BETAIS=0.5;              //two-stream parameter betad for snow
		public static 	double[] EG=new double[]{0.97,0.98};               //emissivity

		
		public RAD_PARAMS()
		{
			
		}
	}
}
