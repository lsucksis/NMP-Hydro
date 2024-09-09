/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2018/1/15
 * Time: 8:38
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;
using System.IO;
using NetCDFDotnet;
using System.Collections.Generic;
using System.Threading.Tasks;
using System.Text;
using System.Threading;


namespace NoahMP
{
	/// <summary>
	/// Description of Driver.
	/// </summary>
	public class Driver
	{
		
		public static string forcePath;
		public static string outputPath;
		public static DateTime startTime = new DateTime(2000, 1, 1, 0, 0, 0);
		public static DateTime endTime;
		/// <summary>
		/// 时间步长：hours
		/// </summary>
		public static TimeSpan timeStep = new TimeSpan(3, 0, 0);
		/// <summary>
		/// 设定1000万步相当于是无限步长
		/// </summary>
		public static int NStep = 10000000;
		/// <summary>
		/// 行数,读入wrfinput文件后，根据HGT的行列数确定
		/// </summary>
		public static int NRows = 0;
		/// <summary>
		/// 列数,读入wrfinput文件后，根据HGT的行列数确定
		/// </summary>
		public static int NCols = 0;
		//public static int NSoil;
		public static string LUTYPE;
		public static int LUCATS = 0;
		public GridCell[,] Cells = null;
		
		/// <summary>
		/// T2D,原变量名T_PHY
		/// </summary>
		public double[,] T3D = null;
		/// <summary>
		/// Q2D
		/// </summary>
		public double[,] QV3D = null;
		/// <summary>
		/// U2D
		/// </summary>
		public double[,] U_PHY = null;
		/// <summary>
		/// V2D
		/// </summary>
		public double[,] V_PHY = null;
		/// <summary>
		/// SWDOWN
		/// </summary>
		public double[,] SWDOWN = null;
		/// <summary>
		/// LWDOWN
		/// </summary>
		public double[,] GLW = null;
		/// <summary>
		/// PSFC
		/// </summary>
		public double[,] P8W3D = null;
		/// <summary>
		/// RAINRATE
		/// </summary>
		public double[,] RAINBL = null;
		public double[,] RUNCOEF = null;
		public double[,] BIAS = null;
		
		public static VEG_PARAMS vegparams = new VEG_PARAMS();
		public void InitModelCells()
		{
			Cells = new GridCell[NRows, NCols];
			for (int row = 0; row < NRows; row++) {
				for (int col = 0; col < NCols; col++) {
					Cells[row, col] = new GridCell();
				}
			}
			
			HFX = new double[NRows, NCols];
			GRDFLX = new double[NRows, NCols];
			ALBEDO = new double[NRows, NCols];
			SNOWC = new double[NRows, NCols];
			CANWAT = new double[NRows, NCols];
			EMISS = new double[NRows, NCols];
			FSAXY = new double[NRows, NCols];
			FIRAXY = new double[NRows, NCols];
			ECANXY = new double[NRows, NCols];
			ETRANXY = new double[NRows, NCols];
			EDIRXY = new double[NRows, NCols];
			
			T2MVXY = new double[NRows, NCols];
			T2MBXY = new double[NRows, NCols];
			Q2MVXY = new double[NRows, NCols];
			// specific humidity to mixing ratio
			Q2MBXY = new double[NRows, NCols];
			// consistent with registry def of Q2
			TRADXY = new double[NRows, NCols];
			NEEXY = new double[NRows, NCols];
			GPPXY = new double[NRows, NCols];
			NPPXY = new double[NRows, NCols];
			FVEGXY = new double[NRows, NCols];
			RUNSFXY = new double[NRows, NCols];
			RUNSBXY = new double[NRows, NCols];
			//double[,]			RUNSBXY =new double[NRows,NCols];
			APARXY = new double[NRows, NCols];
			PSNXY = new double[NRows, NCols];
			SAVXY = new double[NRows, NCols];
			SAGXY = new double[NRows, NCols];
			RSSUNXY = new double[NRows, NCols];
			RSSHAXY = new double[NRows, NCols];
			BGAPXY = new double[NRows, NCols];
			WGAPXY = new double[NRows, NCols];
			TGVXY = new double[NRows, NCols];
			TGBXY = new double[NRows, NCols];
			CHVXY = new double[NRows, NCols];
			CHBXY = new double[NRows, NCols];
			IRCXY = new double[NRows, NCols];
			IRGXY = new double[NRows, NCols];
			SHCXY = new double[NRows, NCols];
			SHGXY = new double[NRows, NCols];
			EVGXY = new double[NRows, NCols];
			GHVXY = new double[NRows, NCols];
			IRBXY = new double[NRows, NCols];
			SHBXY = new double[NRows, NCols];
			EVBXY = new double[NRows, NCols];
			GHBXY = new double[NRows, NCols];
			TRXY = new double[NRows, NCols];
			EVCXY = new double[NRows, NCols];
			CHLEAFXY = new double[NRows, NCols];
			CHUCXY = new double[NRows, NCols];
			CHV2XY = new double[NRows, NCols];
			CHB2XY = new double[NRows, NCols];
			ZLVLXY = new double[NRows, NCols];
			RUNCOEF = new double[NRows, NCols];
			BIAS = new double[NRows, NCols];
			LH = new double[NRows, NCols];
			this.TOTALRUNOFF = new double[NRows, NCols];
			//this.UGDRUNOFF = new double[nrows, ncols];
			this.ACCSFCRUNOFF = new double[NRows, NCols];
			this.ACCUGDRUNOFF = new double[NRows, NCols];
		}
		public Driver()
		{
			//用于确定网格尺寸
			
			
		}
		
		public static double CALC_DECLIN(DateTime NOWDATE, double LATITUDE, double LONGITUDE)
		{
			double DEGRAD = Math.PI / 180;
			double DPD = 360.0 / 365;

			
			//   REAL                           :: HRANG
			//   REAL                           :: DECLIN
			//   REAL                           :: OBECL
			//   REAL                           :: SINOB
			//   REAL                           :: SXLONG
			//   REAL                           :: ARG
			//   REAL                           :: TLOCTIM
			//   INTEGER                        :: IDAY
			//   INTEGER                        :: IHOUR
			//   INTEGER                        :: IMINUTE
			//   INTEGER                        :: ISECOND
			double JULIAN = NOWDATE.DayOfYear - 1 + NOWDATE.Hour / 24.0;
			// CALL GETH_IDTS(NOWDATE(1:10), NOWDATE(1:4)//"-01-01", IDAY)
			//   READ(NOWDATE(12:13), *) IHOUR
			//   READ(NOWDATE(15:16), *) IMINUTE
			//   READ(NOWDATE(18:19), *) ISECOND
			//   JULIAN = REAL(IDAY) + REAL(IHOUR)/24.
			// FOR SHORT WAVE RADIATION

			double DECLIN = 0;

			//-----OBECL : OBLIQUITY = 23.5 DEGREE.

			double OBECL = 23.5 * DEGRAD;
			double SINOB = Math.Sin(OBECL);

			//-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:

			double SXLONG = 0;
			if (JULIAN >= 80) {
				SXLONG = DPD * (JULIAN - 80) * DEGRAD;
			} else {
				SXLONG = DPD * (JULIAN + 285) * DEGRAD;
			}
			double ARG = SINOB * Math.Sin(SXLONG);
			DECLIN = Math.Asin(ARG);

			double TLOCTIM = NOWDATE.Hour + NOWDATE.Minute / 60.0 + NOWDATE.Second / 3600.0 + LONGITUDE / 15.0; // LOCAL TIME IN HOURS
			TLOCTIM = (TLOCTIM + 24.0) % 24;
			double HRANG = 15.0 * (TLOCTIM - 12.0) * DEGRAD;
			double COSZ = Math.Sin(LATITUDE * DEGRAD) * Math.Sin(DECLIN) + Math.Cos(LATITUDE * DEGRAD) * Math.Cos(DECLIN) * Math.Cos(HRANG);
			return COSZ;
		}
		public void ReadParams()
		{
			//该表的内容实际上暂时NoahMP不用它
			StreamReader sr = new StreamReader("VEGPARM.TBL");
			sr.ReadLine();
			LUTYPE = sr.ReadLine();
			string[] line = sr.ReadLine().Split(new char[]{ '\t', ' ', ',' }, StringSplitOptions.RemoveEmptyEntries);
			LUCATS = Convert.ToInt32(line[0]);
			VEG_PARAMS.MVT = LUCATS;
			REDPRM.InitData(LUCATS);
			for (int i = 0; i < LUCATS; i++) {
				line = sr.ReadLine().Split(new char[]{ '\t', ' ', ',' }, StringSplitOptions.RemoveEmptyEntries);
				REDPRM.SHDTBL[i] = Convert.ToDouble(line[1]);
				REDPRM.NROTBL[i] = Convert.ToInt32(line[2]);
				REDPRM.RSTBL[i] = Convert.ToDouble(line[3]);
				REDPRM.RGLTBL[i] = Convert.ToDouble(line[4]);
				REDPRM.HSTBL[i] = Convert.ToDouble(line[5]);
				REDPRM.SNUPTBL[i] = Convert.ToDouble(line[6]);
				REDPRM.MAXALB[i] = Convert.ToDouble(line[7]);
				REDPRM.LAIMINTBL[i] = Convert.ToDouble(line[8]);
				REDPRM.LAIMAXTBL[i] = Convert.ToDouble(line[9]);
				REDPRM.EMISSMINTBL[i] = Convert.ToDouble(line[10]);
				REDPRM.EMISSMAXTBL[i] = Convert.ToDouble(line[11]);
				REDPRM.ALBEDOMINTBL[i] = Convert.ToDouble(line[12]);
				REDPRM.ALBEDOMAXTBL[i] = Convert.ToDouble(line[13]);
				REDPRM.Z0MINTBL[i] = Convert.ToDouble(line[14]);
				REDPRM.Z0MAXTBL[i] = Convert.ToDouble(line[15]);
				REDPRM.ZTOPVTBL[i] = Convert.ToDouble(line[16]);
				REDPRM.ZBOTVTBL[i] = Convert.ToDouble(line[17]);
			}
			sr.Close();
			sr = new StreamReader("SOILPARM.TBL");
			sr.ReadLine();
			string SLTYPE = sr.ReadLine();
			string temp = sr.ReadLine();
			string aaa = temp.Split(new char[]{ ' ', '\t', ',' }, StringSplitOptions.RemoveEmptyEntries)[0];
			NoahMP.NSLTYPE = Convert.ToInt32(aaa);
			for (int i = 0; i < NoahMP.NSLTYPE; i++) {
				string[] values = sr.ReadLine().Split(new char[]{ ' ', '\t', ',' }, StringSplitOptions.RemoveEmptyEntries);
				NoahMP.BB[i] = Convert.ToDouble(values[1]);
				NoahMP.DRYSMC[i] = Convert.ToDouble(values[2]);
				NoahMP.F11[i] = Convert.ToDouble(values[3]);
				NoahMP.MAXSMC[i] = Convert.ToDouble(values[4]);
				NoahMP.REFSMC[i] = Convert.ToDouble(values[5]);
				NoahMP.SATPSI[i] = Convert.ToDouble(values[6]);
				NoahMP.SATDK[i] = Convert.ToDouble(values[7]);
				NoahMP.SATDW[i] = Convert.ToDouble(values[8]);
				NoahMP.WLTSMC[i] = Convert.ToDouble(values[9]);
				NoahMP.QTZ[i] = Convert.ToDouble(values[10]);
			}
			sr.Close();
			sr = new StreamReader("GENPARM.TBL");
			sr.ReadLine();
			string line0 = sr.ReadLine();
			if (line0.Contains("SLOPE_DATA")) {
				REDPRM.SLOPE_DATA = new double[10];
				for (int i = 0; i < 10; i++) {
					REDPRM.SLOPE_DATA[i] = Convert.ToDouble(sr.ReadLine());
				}
			}
			line0 = sr.ReadLine();
			REDPRM.SBETA_DATA = Convert.ToDouble(sr.ReadLine());
			line0 = sr.ReadLine();
			REDPRM.FXEXP_DATA = Convert.ToDouble(sr.ReadLine());
			line0 = sr.ReadLine();
			REDPRM.CSOIL_DATA = Convert.ToDouble(sr.ReadLine());
			line0 = sr.ReadLine();
			REDPRM.SALP_DATA = Convert.ToDouble(sr.ReadLine());
			line0 = sr.ReadLine();
			REDPRM.REFDK_DATA = Convert.ToDouble(sr.ReadLine());
			line0 = sr.ReadLine();
			REDPRM.REFKDT_DATA = Convert.ToDouble(sr.ReadLine());
			line0 = sr.ReadLine();
			REDPRM.FRZK_DATA = Convert.ToDouble(sr.ReadLine());
			line0 = sr.ReadLine();
			REDPRM.ZBOT_DATA = Convert.ToDouble(sr.ReadLine());
			line0 = sr.ReadLine();
			REDPRM.CZIL_DATA = Convert.ToDouble(sr.ReadLine());
			line0 = sr.ReadLine();
			REDPRM.SMLOW_DATA = Convert.ToDouble(sr.ReadLine());
			line0 = sr.ReadLine();
			REDPRM.SMHIGH_DATA = Convert.ToDouble(sr.ReadLine());
			line0 = sr.ReadLine();
			REDPRM.LVCOEF_DATA = Convert.ToDouble(sr.ReadLine());
			sr.Close();
		}
		System.Collections.Generic.SortedSet<string> outputParams = null;
		public void ReadOutputParams()
		{
			outputParams = new SortedSet<string>();
			StreamReader sr = new StreamReader("outputParams.py");
			while (!sr.EndOfStream) {
				string line = sr.ReadLine();
				if (line == "" || line.Trim() == "" || line.Trim().Substring(0, 2) == "//" || line.Trim()[0] == '#')
					continue;
				string param = line.Split(new char[]{ ',' }, StringSplitOptions.RemoveEmptyEntries)[0];
				outputParams.Add(param);
			}
			sr.Close();
		}
		
		public double[] Transform(string[] strs)
		{
			double[] temp = new double[strs.Length];
			for (int i = 0; i < strs.Length - 1; i++) {
				temp[i] = Convert.ToDouble(strs[i + 1]);
			}
			return temp;
			
		}
		public void ReadMPTable(string DatasetName)
		{
			StreamReader sr = new StreamReader("MPTABLE.TBL");
			
			while (true) {
				string line = sr.ReadLine();
				if (line.ToLower().Contains(DatasetName.ToLower()))
					break;
			}
			//sr.ReadLine();
			//string category = sr.ReadLine();
			string description = sr.ReadLine();
			string[] ss = null;
			if (!description.Contains("description".ToUpper())) { //遇到IGBP/MODIS的情形，这样处理
				ss = description.Split(new char[]{ ' ', ',', '\t', '=' }, StringSplitOptions.RemoveEmptyEntries);
			} else { //遇到USGS的情形，这样处理
				ss = sr.ReadLine().Split(new char[]{ ' ', ',', '\t', '=' }, StringSplitOptions.RemoveEmptyEntries);
			}
			VEG_PARAMS.MVT = Convert.ToInt32(ss[1]);
			sr.ReadLine();
			while (true) {
				string line = sr.ReadLine().Trim();
				if (line == "" || line[0] == '!')//(!line.Contains("=")) ||
					continue;
				
				string[] values = line.Split(new char[]{ ' ', ',', '\t', '=' }, StringSplitOptions.RemoveEmptyEntries);
				
				if (values[0] == "ISURBAN")
					vegparams.ISURBAN = Convert.ToInt32(values[1]);
				if (values[0] == "ISWATER")
					vegparams.ISWATER = Convert.ToInt32(values[1]);
				if (values[0] == "ISBARREN")
					vegparams.ISBARREN = Convert.ToInt32(values[1]);
				if (values[0] == "ISSNOW")
					vegparams.ISSNOW = Convert.ToInt32(values[1]);
				if (values[0] == "EBLFOREST")
					vegparams.EBLFOREST = Convert.ToInt32(values[1]);
				if (values[0] == "CH2OP")
					vegparams.CH2OP = Transform(values);
				if (values[0] == "DLEAF")
					vegparams.DLEAF = Transform(values);
				if (values[0] == "Z0MVT")
					vegparams.Z0MVT = Transform(values);
				if (values[0] == "HVT")
					vegparams.HVT = Transform(values);
				if (values[0] == "HVB")
					vegparams.HVB = Transform(values);
				if (values[0] == "DEN")
					vegparams.DEN = Transform(values);
				if (values[0] == "RC")
					vegparams.RC = Transform(values);
				if (values[0] == "RHOL") {
					for (int i = 0; i < values.Length - 1; i++) {
						vegparams.RHOL[i, 0] = Convert.ToDouble(values[i + 1]);
					}
					for (int m = 1; m < 2; m++) {
						line = sr.ReadLine();
						values = line.Split(new char[]{ ' ', ',', '\t', '=' }, StringSplitOptions.RemoveEmptyEntries);
						for (int i = 0; i < values.Length; i++) {
							vegparams.RHOL[i, m] = Convert.ToDouble(values[i]);
						}
					}
				}
				if (values[0] == "RHOS") {
					for (int i = 0; i < values.Length - 1; i++) {
						vegparams.RHOS[i, 0] = Convert.ToDouble(values[i + 1]);
					}
					for (int m = 1; m < 2; m++) {
						line = sr.ReadLine();
						values = line.Split(new char[]{ ' ', ',', '\t', '=' }, StringSplitOptions.RemoveEmptyEntries);
						for (int i = 0; i < values.Length; i++) {
							vegparams.RHOS[i, m] = Convert.ToDouble(values[i]);
						}
					}
				}
				if (values[0] == "TAUL") {
					for (int i = 0; i < values.Length - 1; i++) {
						vegparams.TAUL[i, 0] = Convert.ToDouble(values[i + 1]);
					}
					for (int m = 1; m < 2; m++) {
						line = sr.ReadLine();
						values = line.Split(new char[]{ ' ', ',', '\t', '=' }, StringSplitOptions.RemoveEmptyEntries);
						for (int i = 0; i < values.Length; i++) {
							vegparams.TAUL[i, m] = Convert.ToDouble(values[i]);
						}
					}
				}
				if (values[0] == "TAUS") {
					for (int i = 0; i < values.Length - 1; i++) {
						vegparams.TAUS[i, 0] = Convert.ToDouble(values[i + 1]);
					}
					for (int m = 1; m < 2; m++) {
						line = sr.ReadLine();
						values = line.Split(new char[]{ ' ', ',', '\t', '=' }, StringSplitOptions.RemoveEmptyEntries);
						for (int i = 0; i < values.Length; i++) {
							vegparams.TAUS[i, m] = Convert.ToDouble(values[i]);
						}
					}
				}
				if (values[0] == "XL")
					vegparams.XL = Transform(values);
				if (values[0] == "CWPVT")
					vegparams.CWPVT = Transform(values);
				if (values[0] == "C3PSN")
					vegparams.C3PSN = Transform(values);
				if (values[0] == "KC25")
					vegparams.KC25 = Transform(values);
				if (values[0] == "AKC")
					vegparams.AKC = Transform(values);
				if (values[0] == "KO25")
					vegparams.KO25 = Transform(values);
				if (values[0] == "AKO")
					vegparams.AKO = Transform(values);
				if (values[0] == "AVCMX")
					vegparams.AVCMX = Transform(values);
				if (values[0] == "AQE")
					vegparams.AQE = Transform(values);
				// XL, CWPVT, C3PSN, KC25, AKC, KO25, AKO, AVCMX, AQE,
				if (values[0] == "LTOVRC")
					vegparams.LTOVRC = Transform(values);
				if (values[0] == "DILEFC")
					vegparams.DILEFC = Transform(values);
				if (values[0] == "DILEFW")
					vegparams.DILEFW = Transform(values);
				if (values[0] == "RMF25")
					vegparams.RMF25 = Transform(values);
				if (values[0] == "SLA")
					vegparams.SLA = Transform(values);
				if (values[0] == "FRAGR")
					vegparams.FRAGR = Transform(values);
				if (values[0] == "TMIN")
					vegparams.TMIN = Transform(values);
				if (values[0] == "VCMX25")
					vegparams.VCMX25 = Transform(values);
				if (values[0] == "TDLEF")
					vegparams.TDLEF = Transform(values);
				if (values[0] == "BP")
					vegparams.BP = Transform(values);
				if (values[0] == "MP")
					vegparams.MP = Transform(values);
				if (values[0] == "QE25")
					vegparams.QE25 = Transform(values);
				if (values[0] == "RMS25")
					vegparams.RMS25 = Transform(values);
				if (values[0] == "RMR25")
					vegparams.RMR25 = Transform(values);
				if (values[0] == "ARM")
					vegparams.ARM = Transform(values);
				if (values[0] == "FOLNMX")
					vegparams.FOLNMX = Transform(values);
				if (values[0] == "WDPOOL")
					vegparams.WDPOOL = Transform(values);
				if (values[0] == "WRRAT")
					vegparams.WRRAT = Transform(values);
				if (values[0] == "MRP")
					vegparams.MRP = Transform(values);
				if (values[0] == "SAIM") {
					for (int i = 0; i < values.Length - 1; i++) {
						vegparams.SAIM[i, 0] = Convert.ToDouble(values[i + 1]);
					}
					for (int m = 1; m < 12; m++) {
						line = sr.ReadLine();
						values = line.Split(new char[]{ ' ', ',', '\t', '=' }, StringSplitOptions.RemoveEmptyEntries);
						for (int i = 0; i < values.Length; i++) {
							vegparams.SAIM[i, m] = Convert.ToDouble(values[i]);
						}
					}
				}
				if (values[0] == "LAIM") {
					for (int i = 0; i < values.Length - 1; i++) {
						vegparams.LAIM[i, 0] = Convert.ToDouble(values[i + 1]);
					}
					for (int m = 1; m < 12; m++) {
						line = sr.ReadLine();
						values = line.Split(new char[]{ ' ', ',', '\t', '=' }, StringSplitOptions.RemoveEmptyEntries);
						for (int i = 0; i < values.Length; i++) {
							vegparams.LAIM[i, m] = Convert.ToDouble(values[i]);
						}
					}
				}
				if (values[0] == "SLAREA")
					vegparams.slarea = Transform(values);
				
				if (values[0] == "EPS") {
					vegparams.eps = new double[VEG_PARAMS.MVT, 5];
					for (int i = 0; i < values.Length - 1; i++) {
						vegparams.eps[i, 0] = Convert.ToDouble(values[i + 1]);
					}
					for (int m = 1; m < 5; m++) {
						line = sr.ReadLine();
						values = line.Split(new char[]{ ' ', ',', '\t', '=' }, StringSplitOptions.RemoveEmptyEntries);
						for (int i = 0; i < values.Length; i++) {
							vegparams.eps[i, m] = Convert.ToDouble(values[i]);
						}
					}
					break; //读完EPS就可以结束读取了
				}
				
			}//end while
		}
		public int[,] GetNCVarAsInt(NCFile ncfile, string varName)
		{
//			Console.WriteLine("Reading " + varName);
			int[,] result = null;
			for (int i = 0; i < ncfile.variables.Length; i++) {
				NCVariable var0 = ncfile.variables[i];
				if (var0.name != varName)
					continue;
				int nrow = 0;
				int ncol = 0;
				if (var0.Dimensions.Length == 2) {
					nrow = var0.Dimensions[0].length;
					ncol = var0.Dimensions[1].length;
				} else {
					
					nrow = var0.Dimensions[1].length;
					ncol = var0.Dimensions[2].length;
				}
				result = new int[nrow, ncol];
				if (var0.name == varName) {
					float[] values = var0.ValuesFloat;
					for (int row = 0; row < nrow; row++) {
						for (int col = 0; col < ncol; col++) {
							result[row, col] = (int)values[row * ncol + col];
						}
					}
					return result;
				}
			}
			return result;
		}
		public double[,] GetNCVar(NCFile ncfile, string varName)
		{
			//Console.WriteLine("Reading " + varName);
			double[,] result = null;
			for (int i = 0; i < ncfile.variables.Length; i++) {
				NCVariable var0 = ncfile.variables[i];
				if (var0.name != varName)
					continue;
				int nrow = 0;
				int ncol = 0;
				if (var0.Dimensions.Length == 2) {
					nrow = var0.Dimensions[0].length;
					ncol = var0.Dimensions[1].length;
				} else if (var0.Dimensions.Length == 3) {
					
					nrow = var0.Dimensions[1].length;
					ncol = var0.Dimensions[2].length;
				} else {
					nrow = var0.Dimensions[2].length;
					ncol = var0.Dimensions[3].length;
				}
				result = new double[nrow, ncol];
				if (var0.name == varName) {
					float[] values = var0.ValuesFloat;
					for (int row = 0; row < nrow; row++) {
						for (int col = 0; col < ncol; col++) {
							result[row, col] = (double)values[row * ncol + col];
						}
					}
					return result;
				}
			}
			return result;
		}
		public int[,] GetNCVar_Int(NCFile ncfile, string varName)
		{
			Console.WriteLine("Reading " + varName);
			int[,] result = null;
			for (int i = 0; i < ncfile.variables.Length; i++) {
				NCVariable var0 = ncfile.variables[i];
				if (var0.name != varName)
					continue;
				int ntime = 0;
				int nrow = 0;
				int ncol = 0;
				if (var0.Dimensions.Length == 3) {
					ntime = var0.Dimensions[0].length;
					nrow = var0.Dimensions[1].length;
					ncol = var0.Dimensions[2].length;
				} else if (var0.Dimensions.Length == 2) {
					nrow = var0.Dimensions[0].length;
					ncol = var0.Dimensions[1].length;
				}
				result = new int[nrow, ncol];
				if (var0.name == varName) {
					int[] values = var0.ValuesInt32;
					for (int row = 0; row < nrow; row++) {
						for (int col = 0; col < ncol; col++) {
							result[row, col] = values[row * ncol + col];
						}
					}
					return result;
				}
			}
			return result;
		}
		public double[,,] GetNCVar3D_Restart(NCFile ncfile, string varName)
		{
			Console.WriteLine("Reading " + varName);
			double[,,] result = null;
			for (int i = 0; i < ncfile.variables.Length; i++) {
				NCVariable var0 = ncfile.variables[i];
				if (var0.name != varName)
					continue;
				int ntime = 0;
				int nrow = 0;
				int nlayer = 0;
				int ncol = 0;
				if (var0.Dimensions.Length == 4) {
					ntime = var0.Dimensions[0].length;
					nrow = var0.Dimensions[1].length;
					nlayer = var0.Dimensions[2].length;
					ncol = var0.Dimensions[3].length;
				} else if (var0.Dimensions.Length == 3) {
					nrow = var0.Dimensions[0].length;
					nlayer = var0.Dimensions[1].length;
					ncol = var0.Dimensions[2].length;
				}
				result = new double[nrow, nlayer, ncol];
				if (var0.name == varName) {
					float[] values = var0.ValuesFloat;
					
					for (int row = 0; row < nrow; row++) {
						for (int layer = 0; layer < nlayer; layer++) {
							for (int col = 0; col < ncol; col++) {
								result[row, layer, col] = (double)values[row * ncol * nlayer + layer * ncol + col];
							}
						}
					}
					return result;
				}
			}
			return result;
		}
		public double[,,] GetNCVar3D(NCFile ncfile, string varName)
		{
			Console.WriteLine("Reading " + varName);
			double[,,] result = null;
			for (int i = 0; i < ncfile.variables.Length; i++) {
				NCVariable var0 = ncfile.variables[i];
				if (var0.name != varName)
					continue;
				int time = var0.Dimensions[0].length;
				int nlayer = var0.Dimensions[1].length;
				int nrow = var0.Dimensions[2].length;
				int ncol = var0.Dimensions[3].length;
				result = new double[nrow, nlayer, ncol];
				if (var0.name == varName) {
					float[] values = var0.ValuesFloat;
					
					for (int row = 0; row < nrow; row++) {
						for (int layer = 0; layer < nlayer; layer++) {
							for (int col = 0; col < ncol; col++) {
								double temp = (double)values[nrow * ncol * layer + row * ncol + col];
								result[row, layer, col] =	temp;
							}
						}
					}
					return result;
				}
			}
			return result;
		}
		public void AssignData(double[,] data, string Name)
		{
			
			int nrows = data.GetLength(0);
			int ncols = data.GetLength(1);
			
			switch (Name) {
				case "SOIL_T":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								////Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SNOW_T":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								////Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SMC":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SH2O":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ZSNSO":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SNICE":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SNLIQ":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "QSNOW":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "FWET":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SNEQVO":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "EAH":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "TAH":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ALBOLD":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "CM":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "CH":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ISNOW":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "CANLIQ":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "CANICE":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SNEQV":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SNOWH":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "TV":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "TG":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ZWT":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "WA":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "WT":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "WSLAKE":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "LFMASS":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "RTMASS":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "STMASS":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "WOOD":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "STBLCP":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "FASTCP":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "LAI":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SAI":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "FPAR":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "GVFMIN":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SHDMAX":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "GVFMAX":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ACMELT":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ACSNOW":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "TAUSS":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "QSFC":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SFCRUNOFF":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "UDRUNOFF":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ACCPRCP":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ACCECAN":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ACCEDIR":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ACCETRAN":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SMOISEQ":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "AREAXY":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SMCWTDXY":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "QRFXY":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "DEEPRECHXY":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
					;
			// "SMOISEQ"   , SMOISEQ
			//      "AREAXY"    , AREAXY
			//      "SMCWTDXY"  , SMCWTDXY
			//      "QRFXY"     , QRFXY
			//      "DEEPRECHXY", DEEPRECHXY
			//      "QSPRINGXY" , QSPRINGXY
			//      "QSLATXY"   , QSLATXY
			//      "QRFSXY"    , QRFSXY
			//      "QSPRINGSXY", QSPRINGSXY
			//      "RECHXY"    , RECHXY
			//      "FDEPTHXY"   ,FDEPTHXY
			//      "RIVERCONDXY",RIVERCONDXY
			//      "RIVERBEDXY" ,RIVERBEDXY
			//      "EQZWT"      ,EQZWT
			//      "PEXPXY"     ,PEXPXY
			}
		}
		public void AssignData3D(double[,,] data, string Name)
		{
			
			int nrows = data.GetLength(0);
			int nlayer = data.GetLength(1);
			int ncols = data.GetLength(2);
			
			switch (Name) {
				case "SOIL_T":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								Cells[row, col].STC = new FortDoubleArray(-NSnow + 1, NSoil); //new double[nlayer];
								for (int layer = 0; layer < nlayer; layer++)
									Cells[row, col].STC[layer - NSnow + 1] = data[row, layer, col];
							}
						}
						break;
					}
				case "SNOW_T":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								////Cells[row,col].STC=data[row,col];
								// temperaly do nothing
								//Cells[row,col].sn
							}
						}
						break;
					}
				case "SMC":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SH2O":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ZSNSO":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SNICE":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SNLIQ":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "QSNOW":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "FWET":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SNEQVO":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "EAH":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "TAH":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ALBOLD":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "CM":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "CH":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ISNOW":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "CANLIQ":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "CANICE":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SNEQV":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SNOWH":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "TV":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "TG":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ZWT":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "WA":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "WT":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "WSLAKE":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "LFMASS":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "RTMASS":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "STMASS":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "WOOD":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "STBLCP":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "FASTCP":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "LAI":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SAI":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "FPAR":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "GVFMIN":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SHDMAX":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "GVFMAX":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ACMELT":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ACSNOW":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "TAUSS":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "QSFC":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SFCRUNOFF":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "UDRUNOFF":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ACCPRCP":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ACCECAN":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ACCEDIR":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "ACCETRAN":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SMOISEQ":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "AREAXY":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "SMCWTDXY":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "QRFXY":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
				case "DEEPRECHXY":
					{
						for (int row = 0; row < nrows; row++) {
							for (int col = 0; col < ncols; col++) {
								//Cells[row,col].STC=data[row,col];
							}
						}
						break;
					}
					;
			// "SMOISEQ"   , SMOISEQ
			//      "AREAXY"    , AREAXY
			//      "SMCWTDXY"  , SMCWTDXY
			//      "QRFXY"     , QRFXY
			//      "DEEPRECHXY", DEEPRECHXY
			//      "QSPRINGXY" , QSPRINGXY
			//      "QSLATXY"   , QSLATXY
			//      "QRFSXY"    , QRFSXY
			//      "QSPRINGSXY", QSPRINGSXY
			//      "RECHXY"    , RECHXY
			//      "FDEPTHXY"   ,FDEPTHXY
			//      "RIVERCONDXY",RIVERCONDXY
			//      "RIVERBEDXY" ,RIVERBEDXY
			//      "EQZWT"      ,EQZWT
			//      "PEXPXY"     ,PEXPXY
			}
		}
		
		/*********************************************************************
		 *    reading data from restart file
		 * 
		 * 
		 * 
		 * 
		 *********************************************************************/

		/// <summary>
		/// Soil temperature
		/// </summary>
		public double[,,] TSLB;
		/// <summary>
		/// 雪层温度
		/// </summary>
		public double[,,] TSNOXY;
		public double[,,] SMOIS;
		public double[,,] SH2O;
		public double[,,] ZSNSOXY;
		public double[,,] SNICEXY;
		public double[,,] SNLIQXY;
		public double[,] QSNOWXY;
		public double[,]	FWETXY;
		public double[,]	SNEQVOXY;
		public double[,]	EAHXY;
		public double[,]	TAHXY;
		public double[,]	ALBOLDXY;
		public double[,]	CMXY;
		public double[,]	CHXY;
		/// <summary>
		/// 基于新版本，即0-2表示三层雪。原有的fortran版本默认值为-2,-1,0
		/// </summary>
		public int[,]	ISNOWXY;
		public double[,]	CANLIQXY;
		public double[,]	CANICEXY;
		public double[,]	SNOW;
		public double[,]	SNOWH;
		public double[,]	TVXY;
		public double[,]	TGXY;
		public double[,]	ZWTXY;
		public double[,]	WAXY;
		public double[,]	WTXY;
		public double[,]	WSLAKEXY;
		public double[,]	LFMASSXY;
		public double[,]	RTMASSXY;
		public double[,]	STMASSXY;
		public double[,]	WOODXY;
		public double[,]	STBLCPXY;
		public double[,]	FASTCPXY;
		public double[,]	LAIXY;
		public double[,]	XSAIXY;
		public double[,]	VEGFRA;
		public double[,]	GVFMIN;
		/// <summary>
		/// 等同于是VEGMAX
		/// </summary>
		public double[,]	SHDMAX;
		public double[,] SHDMIN;
		public double[,]	GVFMAX;
		public double[,]	ACSNOM;
		public double[,]	ACSNOW;
		public double[,]	TAUSSXY;
		public double[,] PSFC;
		public double[,]	QSFC;
		/// <summary>
		/// 输出runoff网格，mm in 3h
		/// </summary>
		public double[,]	TOTALRUNOFF;
		
		public double[,]	ACCSFCRUNOFF;
		public double[,] ACCUGDRUNOFF;
		public double[,]	ACCPRCP;
		public double[,]	ACCECAN;
		public double[,]	ACCEDIR;
		public double[,]	ACCETRAN;
		/// <summary>
		/// SMOISEQ[row, z, col]
		/// </summary>
		public double[,,]	SMOISEQ;
		public double[,]	AREAXY;
		public double[,]	SMCWTDXY;
		public double[,]	QRFXY;
		public double[,]	DEEPRECHXY;
		public double[,]	QSPRINGXY;
		public double[,]	QSLATXY;
		public double[,]	QRFSXY;
		public double[,]	QSPRINGSXY;
		public double[,]	RECHXY;
		public double[,]	FDEPTHXY;
		public double[,]	RIVERCONDXY;
		public double[,]	RIVERBEDXY;
		public double[,]	EQZWT;
		public double[,]	PEXPXY;
		
		public double[,] SMSTAV;
		public double[,] SMSTOT;
		
		public double[,] QFX;
		public double[,] LH;
		
		public bool Restart = false;
		
		public static string RestartFile = null;
		/// <summary>
		/// wrfinput文件的路径
		/// </summary>
		public static string LandFile = "wrfinput_d02";
		public static DateTime time0;
		/// <summary>
		/// 0:SingleBox,单点测试模式; 1:Serial multiple box,多点测试模式; 2: Parallel multiple box, 并行模式，有格点输出；3:Parallel multiple box, 并行模式，无格点输出
		/// </summary>
		public static int SingleBoxTest = 2;
		int testRow = 42;
		//93;
		int testCol = 71;
		//184;
		/// <summary>
		/// 汇流模型
		/// </summary>
		RiverRouting routing;
		public void ModifyVegType()
		{
			for (int row = 0; row < NRows; row++) {
				for (int col = 0; col < NCols; col++) {
					if (IVGTYP[row, col] > VEG_PARAMS.MVT)
						IVGTYP[row, col] = vegparams.ISWATER;
				}
			}			
		}
		public void InitStatus()
		{
			Console.WriteLine("Running model start");
			this.ReadNameList();
			this.ReadParams();
			this.ReadOutputParams();
			this.ReadLand(LandFile);
			if (LandFile.Contains("USGS")) {
				this.ReadMPTable("USGS"); //选USGS或者IGBP
			} else if (LandFile.Contains("IGBP")) {
				this.ReadMPTable("IGBP"); //选USGS或者IGBP
				ModifyVegType();
			} else {
				throw new Exception("你的wrfinput文件名中必须包含一个土地利用类型标识：'USGS'或'IGBP'，注意二者的差别。");
			}
			//下一句是刘永和设计用来调参数的，所以要格外小心
			//this.TuneParams();
			this.InitModelCells();
		}
		public void RunModel()
		{
			//初始化汇流模型
			routing = new RiverRouting();
			
			Console.WriteLine("Init the model variables");
			StreamReader sr = new StreamReader("LastRestart.txt");
			string filename0 = sr.ReadLine();
			sr.Close();
			Console.WriteLine("Do you want to run the model based on the last simulation(" + filename0 + ")? (Y or N):");
			string key = Console.ReadLine();
			if (key == "Y" || key == "y") { //以继续上次最后一次模拟
				string year = filename0.Substring(8, 4);
				string month = filename0.Substring(12, 2);
				string day = filename0.Substring(14, 2);
				string hour = filename0.Substring(16, 2);
				RestartFile = filename0;
				time0 = new DateTime(Convert.ToInt32(year), Convert.ToInt32(month), Convert.ToInt32(day), Convert.ToInt32(hour), 0, 0);
				time0 += timeStep;
				this.ReadRestart();
				this.InitRestart();
				
				this.AssignToCell_Restart();
				string routingStatusFile = Path.GetFileNameWithoutExtension(RestartFile) + ".bin";
				if (File.Exists(routingStatusFile)) {
					routing.ReadStatus(routingStatusFile);
				}
			} else {			
				if (Restart) {					//如果是重启续模拟，则需要从文件读入这些状态变量
					this.ReadRestart();
					this.InitRestart();
				
					this.AssignToCell_Restart();
					string routingStatusFile = Path.GetFileNameWithoutExtension(RestartFile) + ".bin";
					if (File.Exists(routingStatusFile)) {
						routing.ReadStatus(routingStatusFile);
					}
					time0 = startTime + timeStep;
				} else {
					this.InitXYParameters();
					this.AssignToCell_NoRestart();
					time0=startTime;
				}
				
			}
			
			for (int row = 0; row < NRows; row++) {
				for (int col = 0; col < NCols; col++) {
					GridCell cell =	Cells[row, col];
					if (cell.ZSNSO == null) {
						Console.WriteLine("No data=" + row + " " + col);
					}
				}
			}
			
			//更新模拟
			//time0 = startTime;
			for (int i = 0; i < NStep; i++) {
				if (time0 >= endTime)
					break;
				Console.WriteLine("Run the step " + i.ToString() + " " + time0);
				string filename = "";
				while (true) {
					filename = time0.Year.ToString("0000") + time0.Month.ToString("00") + time0.Day.ToString("00") + time0.Hour.ToString("00") + ".LDASIN_DOMAIN1";
					if (File.Exists(forcePath + "/" + filename))
						break;
					filename = time0.Year.ToString("0000") + "/" + filename;
					if (File.Exists(forcePath + "/" + filename))
						break;
					
					Console.WriteLine("文件" + forcePath + "/" + filename + "不存在，我先睡一分钟，请您现在准备好该文件。当我发现有文件时我会向后执行。");
					Thread.Sleep(new TimeSpan(0, 1, 0));
					
				}
//				if (time0.Month > 6) {
//					Console.WriteLine("");
//				}
				
				this.ReadForce(forcePath + "/" + filename);
				//陆面模拟
				if (SingleBoxTest < 2) {  //如果是单格子测试，则后续步骤不用进行了
					Console.WriteLine("目前运行在单点调试模式下，请小心使用");
					if (SingleBoxTest == 0)
						Thread.Sleep(100);
					for (int row = 0; row < NRows; row++) {
						//Parallel.For(0, NCols, col =>{
						for (int col = 0; col < NCols; col++) {							
							if (SingleBoxTest == 0 && (row != testRow || col != testCol)) {   //93,184
								continue;
							}
//							Console.WriteLine(row+" "+col);
							GridCell cell =	Cells[row, col];
							if (XLAND[row, col] > 1) {   //对水域不模拟
								ACCSFCRUNOFF[row, col] = -1e36; //下面这些句是自己加的，直至continue
								ACCUGDRUNOFF[row, col] = -1e36;
								LFMASSXY[row, col] = 0;
								RTMASSXY[row, col] = 0;
								STMASSXY[row, col] = 0;
								WOODXY[row, col] = 0;
								continue;
							}
							if (cell.HGT <= 0 || cell.IST == 2)
								continue;
//						if (cell.ICE == 1 || double.IsNaN(cell.TG) || cell.TG < 0
//						    || cell.SH2O[0] > 1e10 || double.IsNaN(cell.SH2O[0])
//						    || cell.NeedCal == false)
//							continue;
//						if(cell.SOILTYP!=1)
//							continue;
							cell.SFCTMP = T3D[row, col];                      // temperature defined at intermediate level [K]
							cell.Q2 = QV3D[row, col] / (1 + QV3D[row, col]);// convert from mixing ratio to specific humidity [kg/kg]
							;
							cell.UU = U_PHY[row, col];                        // u-wind at interface [m/s]
							cell.VV = V_PHY[row, col];                         // v-wind at interface [m/s]
							cell.SOLDN = SWDOWN[row, col];                        // shortwave down from SW scheme [W/m2]
							
							cell.LWDN = GLW[row, col];                           // total longwave down from LW scheme [W/m2]
							cell.SFCPRS = P8W3D[row, col];  // surface pressure defined at intermediate level [Pa]
							//    consistent with temperature, mixing ratio
							cell.PSFC = P8W3D[row, col];                        // surface pressure defined a full levels [Pa]
							cell.PRCP = RAINBL[row, col];                      // timestep precipitation [mm/s]
							//Console.WriteLine("SOLDN "+cell.SOLDN+" SFCTMP "+cell.SFCTMP+" PSFC "+cell.PSFC);
							//下面这句不能加，加上后可能会导致FVEG不及时更新
							//if (VEGFRA != null)
							//	cell.FVEG = VEGFRA[row, col] / 100;
							//Console.WriteLine(row + " " + col);
							if (i == 0 && !Restart) {
								cell.EAH = P8W3D[row, col] * QV3D[row, col] / (0.622 + QV3D[row, col]);
								cell.TAH = T3D[row, col];
								cell.CH = 0.1;
								cell.CM = 0.1;
							}
//						if(i>=121)
//							Console.WriteLine("now");
							this.RunOneStep(row, col, cell);
							
						}
					}					
					time0 += timeStep;
					continue;
				} else {
					for (int row = 0; row < NRows; row++) {
						ParallelOptions option = new ParallelOptions();
						//option.MaxDegreeOfParallelism = 4;
						Parallel.For(0, NCols, col => {
							//for (int col = 0; col < NCols; col++) {
							GridCell cell =	Cells[row, col];
							if (XLAND[row, col] > 1) {   //对水域不模拟
								ACCSFCRUNOFF[row, col] = -1e36; //下面这些句是自己加的，直至continue
								ACCUGDRUNOFF[row, col] = -1e36;
								LFMASSXY[row, col] = 0;
								RTMASSXY[row, col] = 0;
								STMASSXY[row, col] = 0;
								WOODXY[row, col] = 0;
								return;
							}
							if (cell.HGT <= 0 || cell.IST == 2)
								return;
//						if (cell.ICE == 1 || double.IsNaN(cell.TG) || cell.TG < 0
//						    || cell.SH2O[0] > 1e10 || double.IsNaN(cell.SH2O[0])
//						    || cell.NeedCal == false)
//							continue;
//						if(cell.SOILTYP!=1)
//							continue;
							cell.SFCTMP = T3D[row, col];                      // temperature defined at intermediate level [K]
							cell.Q2 = QV3D[row, col] / (1 + QV3D[row, col]);// convert from mixing ratio to specific humidity [kg/kg]
							;
							cell.UU = U_PHY[row, col];                        // u-wind at interface [m/s]
							cell.VV = V_PHY[row, col];                         // v-wind at interface [m/s]
							cell.SOLDN = SWDOWN[row, col];                        // shortwave down from SW scheme [W/m2]
						             	
							cell.LWDN = GLW[row, col];                           // total longwave down from LW scheme [W/m2]
							cell.SFCPRS = P8W3D[row, col];  // surface pressure defined at intermediate level [Pa]
							//    consistent with temperature, mixing ratio
							cell.PSFC = P8W3D[row, col];                        // surface pressure defined a full levels [Pa]
							cell.PRCP = RAINBL[row, col];                      // timestep precipitation [mm/s]
							//Console.WriteLine("SOLDN "+cell.SOLDN+" SFCTMP "+cell.SFCTMP+" PSFC "+cell.PSFC);
//							if (VEGFRA != null)
//								cell.FVEG = VEGFRA[row, col] / 100;
							//Console.WriteLine(row + " " + col);
							if (i == 0 && !Restart) {
								cell.EAH = P8W3D[row, col] * QV3D[row, col] / (0.622 + QV3D[row, col]);
								cell.TAH = T3D[row, col];
								cell.CH = 0.1;
								cell.CM = 0.1;
							}
//						if(i>=121)
//							Console.WriteLine("now");
							this.RunOneStep(row, col, cell);						             	
						});
					}
				}
				this.UpdateMatrix();
				//计算汇流
				if (NoahMP.OPT_UHG != 0)
					routing.RunRouting(time0, TOTALRUNOFF);
				
				TimeSpan span = time0 - startTime;
				if (SingleBoxTest != 3 && span.TotalHours % Output_Frequency_Hours == 0) { //这里刚改过，后续需改为0
					this.WriteOutput();
				}
				//输出重启动文件
				if (span.TotalHours % Restart_Frequency_Hours == 0 && (time0 != startTime)) {
					//this.WriteOutput();
					string filename1 = "RESTART." + time0.Year.ToString("0000") + time0.Month.ToString("00") + time0.Day.ToString("00") + time0.Hour.ToString("00") + "_DOMAIN.nc";
					this.WriteRestart(filename1);
					string filename2 = "RESTART." + time0.Year.ToString("0000") + time0.Month.ToString("00") + time0.Day.ToString("00") + time0.Hour.ToString("00") + "_DOMAIN.bin";
					if (NoahMP.OPT_UHG != 0)
						routing.WriteStatus(filename2);
				}
				time0 += timeStep;
			}
		}
		/// <summary>
		/// 向前运行一步Noah-MP
		/// </summary>
		/// <param name="I">行</param>
		/// <param name="J"></param>
		/// <param name="cell">列</param>
		public void RunOneStep(int I, int J, GridCell cell)
		{
			if (cell.ICE == -1) {   //如果是冰川，则需要按冰川来模拟，但是现在未实现这个

				//     NOAHMP_OPTIONS_GLACIER(IDVEG  ,IOPT_CRS  ,IOPT_BTR  ,IOPT_RUN  ,IOPT_SFC  ,IOPT_FRZ ,
				//    	                       IOPT_INF  ,IOPT_RAD  ,IOPT_ALB  ,IOPT_SNF  ,IOPT_TBOT, IOPT_STC );

				cell.TBOT = Math.Min(cell.TBOT, 263.15);                     // set deep temp to at most -10C
				
				
				//      NOAHMP_GLACIER(       I,       J,    COSZ,   NSNOW,   NSOIL,      DT,  // IN : Time/Space/Model-related
				//                               T_ML,    P_ML,    U_ML,    V_ML,    Q_ML,    SWDN,  // IN : Forcing
				//                               PRCP,    LWDN,    TBOT,    Z_ML, FICEOLD,   ZSOIL,  // IN : Forcing
				//                              QSNOW,  SNEQVO,  ALBOLD,      CM,      CH,   ISNOW,  // IN/OUT :
				//                                SWE,     SMC,   ZSNSO,  SNDPTH,   SNICE,   SNLIQ,  // IN/OUT :
				//                                 TG,     STC,   SMH2O,   TAUSS,  QSFC1D,           // IN/OUT :
				//                                FSA,     FSR,    FIRA,     FSH,    FGEV,   SSOIL,  // OUT :
				//                               TRAD,   ESOIL,   RUNSF,   RUNSB,     SAG,    SALB,  // OUT :
				//                              QSNBOT,PONDING,PONDING1,PONDING2,    T2MB,    Q2MB,  // OUT :
				//                              EMISSI,  FPICE,    CHB2                              // OUT :
				//#ifdef WRF_HYDRO
				//                              , sfcheadrt[I,J]
				//#}
				//);

//							cell.FSNO = 1.0;
//				cell.TV = -1e36;     // Output from standard Noah-MP undefined for glacier points
//				//cell.TGB = TG;
//				cell.CANICE = 0.0;
//				cell.CANLIQ = 0.0;
//				cell.EAH = -1e36;
//				cell.TAH = -1e36;
//				cell.FWET = 0.0;
//				cell.WSLAKE = 0.0;
//				cell.ZWT = -1e36;
//				cell.WA = -1e36;
//				cell.WT = -1e36;
//				cell.LFMASS = 0.0;
//				cell.RTMASS = 0.0;
//				cell.STMASS = 0.0;
//				cell.WOOD = 0.0;
//				cell.STBLCP = -1e36;
//				cell.FASTCP = -1e36;
//				cell.LAI = 0.0;
//				cell.SAI = 0.0;
//				cell.T2MV = -1e36;
//				//cell.Q2MV = -1e36;
//				cell.NEE = 0.0;
//				cell.GPP = 0.0;
//				cell.NPP = 0.0;
//				//cell.FVEGMP = 0.0;
//				cell.ECAN = 0.0;
//				cell.ETRAN = 0.0;
//				cell.APAR = 0.0;
//				cell.PSN = 0.0;
//				cell.SAV = 0.0;
//				cell.RSSUN = -1e36;
//				cell.RSSHA = -1e36;
//				cell.BGAP = -1e36;
//				cell.WGAP = -1e36;
//				cell.TGV = -1e36;
//				cell.CHV = -1e36;
//				//cell.CHB = CH;
//				cell.IRC = -1e36;
//				cell.IRG = -1e36;
//				cell.SHC = -1e36;
//				cell.SHG = -1e36;
//				cell.EVG = -1e36;
//
//				cell.GHV = -1e36;
				////							cell.IRB = FIRA;
				////							cell.SHB = FSH;
				////							cell.EVB = FGEV;
				////							cell.GHB = SSOIL;
//				cell.TR = 0.0;
//				cell.EVC = 0.0;
//				cell.CHLEAF = -1e36;
//				cell.CHUC = -1e36;
//				cell.CHV2 = -1e36;
//				cell.FCEV = 0.0;
//				cell.FCTR = 0.0;
//							cell.QFX[I, J] = ESOIL;
//							cell.LH[I, J] = FGEV;


			} else {
				cell.ILOC = I;
				cell.JLOC = J;
				cell.SurfaceFlux();
				
//				cell.CheckModel();
				LH[I, J] = cell.LH;
				//QFX[I,J] = cell.QFX;
				
				TOTALRUNOFF[I, J] = cell.Q_HillSlope * DT;
				ACCSFCRUNOFF[I, J] += cell.RUNSRF * DT;
				ACCUGDRUNOFF[I, J] += cell.RUNSUB * DT;
//				if(ACCUGDRUNOFF[I, J]>1500)
//					throw new Exception("");
				
				ACCPRCP[I, J] += cell.PRCP * DT;
				ACCECAN[I, J] += cell.ECAN * DT;
				ACCETRAN[I, J] += cell.ETRAN * DT;
				ACCEDIR[I, J] += cell.EDIR * DT;
				ACSNOW[I, J] += cell.PRCP * cell.FPICE;
				ACSNOM[I, J] += cell.QSNBOT * DT + cell.PONDING + cell.PONDING1 + cell.PONDING2;
				cell.sumPcp += cell.PRCP * DT;
				cell.sumEvap += (cell.ECAN + cell.ETRAN + cell.EDIR) * DT;
				cell.sumRunoff += (cell.RUNSRF + cell.RUNSUB) * DT;
				RUNCOEF[I, J] = cell.sumRunoff / cell.sumPcp;
				BIAS[I, J] = (cell.sumEvap + cell.sumRunoff - cell.sumPcp) / cell.sumPcp;
//				if (I == testRow && J == testCol) {   //93,184
//					Console.WriteLine("PRCP=" + cell.sumPcp.ToString("0.0000000000") + " EVAP=" + cell.sumEvap.ToString("0.0000000000") + " Runoff=" + (cell.sumRunoff).ToString("0.0000000000") +
//					                  " Bias=" + ((cell.sumRunoff + cell.sumEvap - cell.sumPcp) / cell.sumPcp).ToString("0.0000") + " Runoff Coef=" + (cell.sumRunoff / cell.sumPcp).ToString("0.0000"));
//				}
			} // glacial split ends
		}
		
		
		
		/*********************************************************************
		 *     Reading data from land file
		 * 
		 *********************************************************************/
		
		public double[,] latitude;
		public double[,] longitude;
		//public static DateTime NowDate;
		/// <summary>
		/// 水域=2，陆地=1，其它我尚未获得相关说明
		/// </summary>
		public double[,] XLAND;
		public double[,] SEAICE;
		public double[,] TERRAIN;
		public double[,] TMN;
		public double[,] MSFTX;
		public double[,] MSFTY;
		public int[,] IVGTYP;
		public int[,] ISLTYP;
		
		public double[,] XICE;
		public double XICE_THRES = 0.5;
		
		public int ISICE = 24;
		public int ISURBAN = 1;
		public int	ISOILWATER = 14;
		public int ISWATER = 16;
		public int ISLAKE = -1;
		public static double DX = 6000;
		public static double DY = 6000;
		/// <summary>
		/// 单位s
		/// </summary>
		public static double DT = 10800;
		public string MAP_PROJ_CHAR = "Lambert Conformal";
		
		
		public void ReadForce(string filename)
		{
			NCFile file = new NCFile(filename, OpenMode.read);
			this.QV3D = GetNCVar(file, "Q2D");
			this.T3D = GetNCVar(file, "T2D");
			this.U_PHY = GetNCVar(file, "U2D");
			this.V_PHY = GetNCVar(file, "V2D");
			this.P8W3D = GetNCVar(file, "PSFC");
			this.GLW = GetNCVar(file, "LWDOWN");
			this.SWDOWN = GetNCVar(file, "SWDOWN");
			this.RAINBL = GetNCVar(file, "RAINRATE");
			double[,] vegf = GetNCVar(file, "VEGFRA");
			if (vegf != null)
				this.VEGFRA = vegf;
			file.Close();
		}
		public void ReadLand(string filename)
		{
			NCFile file = new NCFile(filename, OpenMode.read);
			latitude = GetNCVar(file, "XLAT");
			longitude = GetNCVar(file, "XLONG");
			XLAND = GetNCVar(file, "XLAND");
			SEAICE = GetNCVar(file, "SEAICE");
			TERRAIN = GetNCVar(file, "HGT");
			TMN = GetNCVar(file, "TMN");
			MSFTX = GetNCVar(file, "MAPFAC_MX");
			MSFTY = GetNCVar(file, "MAPFAC_MY");
			try {
				IVGTYP = GetNCVar_Int(file, "IVGTYP");
				ISLTYP = GetNCVar_Int(file, "ISLTYP");
			} catch {
				IVGTYP = GetNCVarAsInt(file, "IVGTYP");
				ISLTYP = GetNCVarAsInt(file, "ISLTYP");
			}
			
			VEGFRA = GetNCVar(file, "VEGFRA");
			LAIXY = GetNCVar(file, "LAI");
			GVFMIN = GetNCVar(file, "GVFMIN");
			SHDMAX = GetNCVar(file, "SHDMAX");
			SHDMIN = GetNCVar(file, "SHDMIN");
			SNOW = GetNCVar(file, "SNOW");
			SNOWH = GetNCVar(file, "SNOWH");
			SMOIS = GetNCVar3D(file, "SMOIS");
			QSFC = GetNCVar(file, "Q2");
			PSFC = GetNCVar(file, "PSFC");
			TSK = GetNCVar(file, "TSK");
			SH2O = GetNCVar3D(file, "SH2O");
			TSLB = GetNCVar3D(file, "TSLB");
			file.Close();
			NRows = TERRAIN.GetLength(0);
			NCols = TERRAIN.GetLength(1);
			Console.WriteLine("The domain size is " + NRows + " rows and " + NCols + " cols.");
			if (XLAND == null) {
				XLAND = new double[NRows, NCols];
			}
			if (MSFTX == null) {
				MSFTX = new double[NRows, NCols];
			}
			if (MSFTY == null) {
				MSFTY = new double[NRows, NCols];
			}
			if (LAIXY == null) {
				LAIXY = new double[NRows, NCols];
			}
			if (GVFMIN == null) {
				GVFMIN = new double[NRows, NCols];
			}
			if (SHDMAX == null) {
				SHDMAX = new double[NRows, NCols];
			}
			if (SHDMIN == null) {
				SHDMIN = new double[NRows, NCols];
			}
			if (SNOW == null) {
				SNOW = new double[NRows, NCols];
			}
			if (QSFC == null) {
				QSFC = new double[NRows, NCols];
			}
			if (PSFC == null) {
				PSFC = new double[NRows, NCols];
			}
		}
		/// <summary>
		/// 与namelist中一致
		/// </summary>
		public static int NSoil = 4;
		public static int NSnow = 3;
		/// <summary>
		/// reference height (m)
		/// </summary>
		public static double ZLVL = 15;
		public static FortDoubleArray DZS = new FortDoubleArray(new double[] {
			0.1,
			0.3,
			0.6,
			1
		});
		public static int KDAY = -1;
		public static int KHOUR = -1;
		public static int Restart_Frequency_Hours = 24;
		public static int Output_Frequency_Hours = 24;
		
		/// <summary>
		/// 读入namelist
		/// </summary>
		public void ReadNameList()
		{
			StreamReader sr = new StreamReader("namelist.hrldas", Encoding.UTF8);
			int start_year = -1;
			int start_month = -1;
			int start_day = -1;
			int start_hour = -1;
//			int KDAY = -1;
//			int KHOUR = -1;
//
			while (!sr.EndOfStream) {
				string s = sr.ReadLine().Trim();
				if (s == "" || s[0] == '!' || s[0] == '&' || s[0] == '/')
					continue;
				string[] line = s.Split(new char[]{ '=' });
				line[0] = line[0].Trim();
				line[1] = line[1].Trim();
				if (line[1].Contains("!")) {
					int pos = line[1].IndexOf("!");
					if (pos > 0)
						line[1] = line[1].Substring(0, pos);
				}
				if (line[0] == "HRLDAS_CONSTANTS_FILE")
					LandFile = line[1];
				if (line[0] == "INDIR") {
					forcePath = line[1];
				}
				if (line[0] == "OUTDIR") {
					outputPath = line[1];
				}
				if (line[0] == "START_YEAR") {
					start_year = Convert.ToInt32(line[1]);
				}
				if (line[0] == "START_MONTH") {
					start_month = Convert.ToInt32(line[1]);
				}
				if (line[0] == "START_DAY") {
					start_day = Convert.ToInt32(line[1]);
				}
				if (line[0] == "START_HOUR") {
					start_hour = Convert.ToInt32(line[1]);
				}
				if (line[0] == "KHOUR") {
					KHOUR = Convert.ToInt32(line[1]);
				}
				if (line[0] == "KDAY") {
					KDAY = Convert.ToInt32(line[1]);
				}
				if (line[0] == "RESTART_FILENAME_REQUESTED") {
					Restart = true;
					RestartFile = line[1];
				}
				if (line[0] == "DYNAMIC_VEG_OPTION") {
					NoahMP.DVEG = Convert.ToInt32(line[1]);
				}
				if (line[0] == "CANOPY_STOMATAL_RESISTANCE_OPTION") {
					NoahMP.OPT_CRS = Convert.ToInt32(line[1]);
				}
				if (line[0] == "BTR_OPTION") {
					NoahMP.OPT_BTR = Convert.ToInt32(line[1]);
				}
				if (line[0] == "RUNOFF_OPTION") {
					NoahMP.OPT_RUN = Convert.ToInt32(line[1]);
				}
				if (line[0] == "SURFACE_DRAG_OPTION") {
					NoahMP.OPT_SFC = Convert.ToInt32(line[1]);
				}
				if (line[0] == "FROZEN_SOIL_OPTION") {
					NoahMP.OPT_INF = Convert.ToInt32(line[1]);
				}
				if (line[0] == "SUPERCOOLED_WATER_OPTION") {
					NoahMP.OPT_FRZ = Convert.ToInt32(line[1]);
				}
				if (line[0] == "RADIATIVE_TRANSFER_OPTION") {
					NoahMP.OPT_RAD = Convert.ToInt32(line[1]);
				}
				if (line[0] == "SNOW_ALBEDO_OPTION") {
					NoahMP.OPT_ALB = Convert.ToInt32(line[1]);
				}
				if (line[0] == "PCP_PARTITION_OPTION") {
					NoahMP.OPT_SNF = Convert.ToInt32(line[1]);
				}
				if (line[0] == "TBOT_OPTION") {
					NoahMP.OPT_TBOT = Convert.ToInt32(line[1]);
				}
				if (line[0] == "TEMP_TIME_SCHEME_OPTION") {
					NoahMP.OPT_STC = Convert.ToInt32(line[1]);
				}
				if (line[0] == "FORCING_TIMESTEP") {
					int value = Convert.ToInt32(line[1]);
					Driver.timeStep = new TimeSpan(0, 0, value);
				}
				if (line[0] == "NOAH_TIMESTEP") {
					int value = Convert.ToInt32(line[1]);
					TimeSpan span = new TimeSpan(0, 0, value);
					//not implemented now
				}
				if (line[0] == "OUTPUT_TIMESTEP") {
					int value = Convert.ToInt32(line[1]);
					Output_Frequency_Hours = value / 3600;
					//not implemented now
				}
				if (line[0] == "RESTART_FREQUENCY_HOURS") {
					Restart_Frequency_Hours = Convert.ToInt32(line[1]);
					//TimeSpan span = new TimeSpan(Restart_Frequency_Hours, 0, 0);
					//not implemented now
				}
				if (line[0] == "NSOIL") {
					NoahMP.NSOIL = Convert.ToInt32(line[1]);
				}
				if (line[0] == "soil_thick_input[1]") {
					Driver.DZS[1] = Convert.ToDouble(line[1]);
				}
				if (line[0] == "soil_thick_input(2)") {
					Driver.DZS[2] = Convert.ToDouble(line[1]);
				}
				if (line[0] == "soil_thick_input(3)") {
					Driver.DZS[3] = Convert.ToDouble(line[1]);
				}
				if (line[0] == "soil_thick_input(4)") {
					Driver.DZS[4] = Convert.ToDouble(line[1]);
				}
				if (line[0] == "ZLVL") {
					Driver.ZLVL = Convert.ToDouble(line[1]);
				}
				if (line[0] == "HILL_SLOPE_RUNOFF_OPTION") {
					NoahMP.OPT_UHG = Convert.ToInt32(line[1]);
				}
				if (line[0] == "threads") {
					NoahMP.threads = Convert.ToInt32(line[1]);
				}
			}
			sr.Close();
			startTime = new DateTime(start_year, start_month, start_day, start_hour, 0, 0);
			//endTime = new DateTime(end_year, end_month, end_day, end_hour, 0, 0);
			if (KHOUR > 0) {
				TimeSpan span = new TimeSpan(KHOUR, 0, 0);
				endTime = startTime + span;
			} else {
				if (KDAY < 0)
					throw new Exception("Time span is set incorrectly");
				TimeSpan span = new TimeSpan(KDAY, 0, 0, 0);
				endTime = startTime + span;
			}
		}
		public void ReadInput(string filename)
		{
//			SortedList<string,string> attrs = new SortedList<string, string>();
//			NCFile file = new NCFile(filename, OpenMode.read);
//			for (int i = 0; i < file.attributes.Length; i++) {
//				string name = file.attributes[i].name;
//				attrs[name] =	file.attributes[i].value;
//			}
			
//			NSOIL=Convert.ToInt32(attrs["soil_layers_stag"]);
//			NSNOW=Convert.ToInt32(attrs["snow_layers"]);
//			NCols=Convert.ToInt32(attrs["west_east"]);
//			NRows=Convert.ToInt32(attrs["south_north"]);

		}
		public void SnowInit(int NSNOW, int NSOIL, double[] ZSOIL, double[,] SWE, double[,] TGXY, double[,] SNOWH,
			double[,,] ZSNSOXY, double[,,] TSNOXY, double[,,] SNICEXY, double[,,] SNLIQXY, int[,] ISNOWXY)
		{
			double[] DZSNO = new double[NSnow];
			double[] DZSNSO = new double[NSnow + NSoil];
			for (int J = 0; J < NCols; J++) {
				for (int I = 0; I < NRows; I++) {
					if (I != 23 || J != 38)
						continue;
					if (SNOWH[I, J] < 0.025) {
						ISNOWXY[I, J] = 0;
						for (int i = 0; i < NSnow; i++) {
							DZSNO[i] = 0;
						}
					} else {
						
						if ((SNOWH[I, J] >= 0.025) && (SNOWH[I, J] <= 0.05)) {
							ISNOWXY[I, J] = NSNOW - 2;
							DZSNO[NSNOW - 1] = SNOWH[I, J];
						} else if ((SNOWH[I, J] > 0.05) && (SNOWH[I, J] <= 0.10)) {
							ISNOWXY[I, J] = -2;
							DZSNO[NSNOW - 2] = SNOWH[I, J] / 2;
							DZSNO[NSNOW - 1] = SNOWH[I, J] / 2;
						} else if ((SNOWH[I, J] > 0.10) && (SNOWH[I, J] <= 0.25)) {
							ISNOWXY[I, J] = -2;
							DZSNO[NSNOW - 2] = 0.05;
							DZSNO[NSNOW - 1] = SNOWH[I, J] - DZSNO[NSNOW - 2];
						} else if ((SNOWH[I, J] > 0.25) && (SNOWH[I, J] <= 0.45)) {
							ISNOWXY[I, J] = -3;
							DZSNO[NSNOW - 3] = 0.05;
							DZSNO[NSNOW - 2] = 0.5 * (SNOWH[I, J] - DZSNO[NSNOW - 3]);
							DZSNO[NSNOW - 1] = 0.5 * (SNOWH[I, J] - DZSNO[NSNOW - 3]);
						} else if (SNOWH[I, J] > 0.45) {
							ISNOWXY[I, J] = -3;
							DZSNO[NSNOW - 3] = 0.05;
							DZSNO[NSNOW - 2] = 0.20;
							DZSNO[NSNOW - 1] = SNOWH[I, J] - DZSNO[NSNOW - 2] - DZSNO[NSNOW - 3];
						} else
							throw new Exception();
						
						//                write(6,*) "SNODEP[I,J] = ", SNODEP[I,J]
						//                CALL wrf_error_fatal("Problem with the logic assigning snow layers.")

						for (int i = 0; i < NSnow; i++) {
							TSNOXY[I, i, J] = 0;
							SNICEXY[I, i, J] = 0;
							SNLIQXY[I, i, J] = 0;
						}
						int ISNOW = ISNOWXY[I, J];
						for (int IZ = ISNOW; IZ < NSnow; IZ++) {
							
							TSNOXY[I, IZ, J] = TGXY[I, J]; // [k]
							SNLIQXY[I, IZ, J] = 0.00;
							SNICEXY[I, IZ, J] = 1.00 * DZSNO[IZ] * (SNOW[I, J] / SNOWH[I, J]);  // [kg/m3]
						}

						// Assign local variable DZSNSO, the soil/snow layer thicknesses, for snow layers
						for (int IZ = ISNOW; IZ < NSnow; IZ++) {
							DZSNSO[IZ] = -DZSNO[IZ];
						}

						// Assign local variable DZSNSO, the soil/snow layer thicknesses, for soil layers
						DZSNSO[NSNOW] = ZSOIL[0];
						for (int IZ = 1; IZ < NSoil; IZ++) {
							
							DZSNSO[IZ + NSnow] = (ZSOIL[IZ] - ZSOIL[IZ - 1]);
						}

						// Assign ZSNSOXY, the layer depths, for soil and snow layers
						ZSNSOXY[I, ISNOW, J] = DZSNSO[ISNOW];
						for (int IZ = ISNOW + 1; IZ < NSoil; IZ++) {
							ZSNSOXY[I, IZ, J] = ZSNSOXY[I, IZ - 1, J] + DZSNSO[IZ];
						}
					}
				}
				
			}
		}
		public void InitXYParameters()
		{
			double BLIM = 5.5;
			double HLICE = 3.335E5;
			double GRAV = 9.81;
			double T0 = 273.15;
			double XLAND = 1.0;//   // water = 2.0, land = 1.0
			double[,] XICE = new double[NRows, NCols];  //   // fraction of grid that is seaice
			double XICE_THRESHOLD = 0.5;//   ! fraction of grid determining seaice
			
			this.SNOWH = new double[NRows, NCols];
			//this.SNOW.CopyTo(this.SNOWH, 0);
			for (int row = 0; row < NRows; row++) {
				for (int col = 0; col < NCols; col++) {
					this.SNOWH[row, col] *= 0.0005;
				}
			}

			ACCECAN = new double[NRows, NCols];
			ACCEDIR = new double[NRows, NCols];
			ACCETRAN = new double[NRows, NCols];
			ACCPRCP = new double[NRows, NCols];
			ACSNOM = new double[NRows, NCols];
			ACSNOW = new double[NRows, NCols];
			TVXY = new double[NRows, NCols];
			TGXY = new double[NRows, NCols];
			ISNOWXY = new int[NRows, NCols];
			CANLIQXY = new double[NRows, NCols];
			CANICEXY = new double[NRows, NCols];
			EAHXY = new double[NRows, NCols];
			TAHXY = new double[NRows, NCols];
			CMXY = new double[NRows, NCols];
			CHXY = new double[NRows, NCols];
			FWETXY = new double[NRows, NCols];
			SNEQVOXY = new double[NRows, NCols];
			ALBOLDXY = new double[NRows, NCols];
			QSNOWXY = new double[NRows, NCols];
			WSLAKEXY = new double[NRows, NCols];
			ZWTXY = new double[NRows, NCols];
			WAXY = new double[NRows, NCols];
			WTXY = new double[NRows, NCols];
			ZSNSOXY = new double[NRows, NSoil + NSnow, NCols];
			TSNOXY = new double[NRows, NSoil + NSnow, NCols];
			SNICEXY = new double[NRows, NSnow, NCols];
			SNLIQXY = new double[NRows, NSnow, NCols];
			STMASSXY = new double[NRows, NCols];
			LFMASSXY = new double[NRows, NCols];
			WOODXY = new double[NRows, NCols];
			STBLCPXY = new double[NRows, NCols];
			FASTCPXY = new double[NRows, NCols];
			LAIXY = new double[NRows, NCols];
			XSAIXY = new double[NRows, NCols];
			TAUSSXY = new double[NRows, NCols];
			RTMASSXY = new double[NRows, NCols];
			WAXY = new double[NRows, NCols];
			WTXY = new double[NRows, NCols];
			CANWAT = new double[NRows, NCols];
			T2MVXY = new double[NRows, NCols];
			T2MBXY = new double[NRows, NCols];
			RUNSFXY = new double[NRows, NCols];
			
			if (NoahMP.OPT_RUN == 5) {
				AREAXY = new double[NRows, NCols];
			}

			//SH2O=GetNCVar3D_Restart(file, "SH2O");
			
			
//			 DO J = 0,jtf
			//          DO I = its,itf
			for (int J = 0; J < NCols; J++) {
				for (int I = 0; I < NRows; I++) {
					TVXY[I, J] = TSK[I, J];
					if (SNOW[I, J] > 0.0 && TSK[I, J] > 273.15)
						TVXY[I, J] = 273.15;
					TGXY[I, J] = TSK[I, J];
					if (SNOW[I, J] > 0.0 && TSK[I, J] > 273.15)
						TGXY[I, J] = 273.15;
					CANWAT[I, J] = 0.0;
					CANLIQXY[I, J] = CANWAT[I, J];
					CANICEXY[I, J] = 0;
					EAHXY[I, J] = 2000;
					TAHXY[I, J] = TSK[I, J];
					if (SNOW[I, J] > 0.0 && TSK[I, J] > 273.15)
						TAHXY[I, J] = 273.15;
					//             tahxy      [I,J] = 287;
					//jref:start
					T2MVXY[I, J] = TSK[I, J];
					if (SNOW[I, J] > 0.0 && TSK[I, J] > 273.15)
						T2MVXY[I, J] = 273.15;
					T2MBXY[I, J] = TSK[I, J];
					if (SNOW[I, J] > 0.0 && TSK[I, J] > 273.15)
						T2MBXY[I, J] = 273.15;
					//chstarxy     [I,J] = 0.1;
					//jref:end

					CMXY[I, J] = 0.0;
					CHXY[I, J] = 0.0;
					FWETXY[I, J] = 0.0;
					SNEQVOXY[I, J] = 0.0;
					ALBOLDXY[I, J] = 0.65;
					QSNOWXY[I, J] = 0.0;
					WSLAKEXY[I, J] = 0.0;

					if (NoahMP.OPT_RUN != 5) {
						WAXY[I, J] = 4900;
						WTXY[I, J] = WAXY[I, J];                                   //???
						ZWTXY[I, J] = (25 + 2.0) - WAXY[I, J] / 1000 / 0.2;            //???
					} else {
						WAXY[I, J] = 0;
						WTXY[I, J] = 0;
						AREAXY[I, J] = (DX * DY) / (MSFTX[I, J] * MSFTY[I, J]);
					}

					LFMASSXY[I, J] = 50;         //
					STMASSXY[I, J] = 50.0;       //
					RTMASSXY[I, J] = 500.0;     //
					WOODXY[I, J] = 500.0;    //
					STBLCPXY[I, J] = 1000.0;    //
					FASTCPXY[I, J] = 1000.0;    //
					XSAIXY[I, J] = 0.1;        //
					
					RUNSFXY[I, J] = -1e36; //刘永和自己加的

				}
			}
			// initialize soil liquid water content SH2O
			for (int J = 0; J < NCols; J++) {
				for (int I = 0; I < NRows; I++) {
//					if (I != 28 || J != 45)
//							continue;
					if (IVGTYP[I, J] == ISICE && XICE[I, J] <= 0.0) {
						for (int NS = 0; NS < NSoil; NS++) {
							
							SMOIS[I, NS, J] = 1.0;                     // glacier starts all frozen
							SH2O[I, NS, J] = 0.0;
							TSLB[I, NS, J] = Math.Min(TSLB[I, NS, J], 263.15); // set glacier temp to at most -10C
						}
						//TMN[I,J] = MIN(TMN[I,J],263.15);         // set deep temp to at most -10C
						SNOW[I, J] = Math.Max(SNOW[I, J], 10.0);        // set SWE to at least 10mm
						SNOWH[I, J] = SNOW[I, J] * 0.01;               // SNOW in mm and SNOWH in m
					} else {
						
						double BX = NoahMP.BB[ISLTYP[I, J] - 1];
						double SMCMAX = NoahMP.MAXSMC[ISLTYP[I, J] - 1];
						for (int NS = 0; NS < NSoil; NS++) {
							if (SMOIS[I, NS, J] > SMCMAX)
								SMOIS[I, NS, J] = SMCMAX;
						}
						double PSISAT = NoahMP.SATPSI[ISLTYP[I, J] - 1];
						if ((BX > 0.0) && (SMCMAX > 0.0) && (PSISAT > 0.0)) {
							for (int NS = 0; NS < NSoil; NS++) {
								double tslb = TSLB[I, NS, J];
								if (tslb < 273.149) {    // Use explicit as initial soil ice
									double FK = Math.Pow(HLICE / (GRAV * (-PSISAT)) * ((tslb - T0) / tslb), -1 / BX) * SMCMAX;
									FK = Math.Max(FK, 0.02);
									SH2O[I, NS, J] = Math.Min(FK, SMOIS[I, NS, J]);
								} else {
									SH2O[I, NS, J] = SMOIS[I, NS, J];
								}
							}
						} else {
							for (int NS = 0; NS < NSoil; NS++) {
								SH2O[I, NS, J] = SMOIS[I, NS, J];
							}
						}
					}
				}
			}

			double[] ZSOIL = new double[NSoil];
			// Given the soil layer thicknesses (in DZS), initialize the soil layer
			// depths from the surface.
			ZSOIL[0] = -DZS.data[0];         // negative
			for (int NS = 1; NS < NSoil; NS++) {
				ZSOIL[NS] = ZSOIL[NS - 1] - DZS[NS];
			}
			SnowInit(NSnow, NSoil, ZSOIL, SNOW, TGXY, SNOWH, ZSNSOXY, TSNOXY, SNICEXY, SNLIQXY, ISNOWXY);
			//SNOW_INIT
			
		}
		public void ReadRestart()
		{
			Console.WriteLine("Reading restart file:" + RestartFile);
			//NetCDFDotnet.NCFile file = new NCFile("RESTART.1991060100_DOMAIN2", OpenMode.read);
			NetCDFDotnet.NCFile file = new NCFile(RestartFile, OpenMode.read);
			this.TSLB = GetNCVar3D_Restart(file, "SOIL_T");
			this.TSNOXY = GetNCVar3D_Restart(file, "SNOW_T");
			this.SMOIS = GetNCVar3D_Restart(file, "SMC");
			this.SH2O = GetNCVar3D_Restart(file, "SH2O");
			this.ZSNSOXY = GetNCVar3D_Restart(file, "ZSNSO");
			this.SNICEXY = GetNCVar3D_Restart(file, "SNICE");
			this.SNLIQXY = GetNCVar3D_Restart(file, "SNLIQ");
			this.QSNOWXY = GetNCVar(file, "QSNOW");
			this.FWETXY = GetNCVar(file, "FWET");
			this.SNEQVOXY = GetNCVar(file, "SNEQVO");
			this.EAHXY = GetNCVar(file, "EAH");
			this.TAHXY = GetNCVar(file, "TAH");
			this.ALBOLDXY = GetNCVar(file, "ALBOLD");
			this.CMXY = GetNCVar(file, "CM");
			this.CHXY = GetNCVar(file, "CH");
			this.ISNOWXY = GetNCVar_Int(file, "ISNOW");
			this.CANLIQXY = GetNCVar(file, "CANLIQ");
			this.CANICEXY = GetNCVar(file, "CANICE");
			this.SNOW = GetNCVar(file, "SNEQV");
			this.SNOWH = GetNCVar(file, "SNOWH");
			this.TVXY = GetNCVar(file, "TV");
			this.TGXY = GetNCVar(file, "TG");
			this.ZWTXY = GetNCVar(file, "ZWT");
			this.WAXY = GetNCVar(file, "WA");
			this.WTXY = GetNCVar(file, "WT");
			this.WSLAKEXY = GetNCVar(file, "WSLAKE");
			this.LFMASSXY = GetNCVar(file, "LFMASS");
			this.RTMASSXY = GetNCVar(file, "RTMASS");
			this.STMASSXY = GetNCVar(file, "STMASS");
			this.WOODXY = GetNCVar(file, "WOOD");
			this.STBLCPXY = GetNCVar(file, "STBLCP");
			this.FASTCPXY = GetNCVar(file, "FASTCP");
			this.LAIXY = GetNCVar(file, "LAI");
			this.XSAIXY = GetNCVar(file, "SAI");
			this.VEGFRA = GetNCVar(file, "FPAR");
			this.GVFMIN = GetNCVar(file, "GVFMIN");
			this.SHDMAX = GetNCVar(file, "SHDMAX");
			this.SHDMIN = GetNCVar(file, "SHDMIN");
			this.GVFMAX = GetNCVar(file, "GVFMAX");
			this.ACSNOM = GetNCVar(file, "ACMELT");
			this.ACSNOW = GetNCVar(file, "ACSNOW");
			this.TAUSSXY = GetNCVar(file, "TAUSS");
			this.QSFC = GetNCVar(file, "QSFC");
			this.ACCSFCRUNOFF = GetNCVar(file, "SFCRUNOFF");
			this.ACCUGDRUNOFF = GetNCVar(file, "UGDRUNOFF");
			this.ACCPRCP = GetNCVar(file, "ACCPRCP");
			this.ACCECAN = GetNCVar(file, "ACCECAN");
			this.ACCEDIR = GetNCVar(file, "ACCEDIR");
			this.ACCETRAN = GetNCVar(file, "ACCETRAN");
			this.SMOISEQ = GetNCVar3D_Restart(file, "SMOISEQ");
			this.AREAXY = GetNCVar(file, "AREAXY");
			this.SMCWTDXY = GetNCVar(file, "SMCWTDXY");
			this.QRFXY = GetNCVar(file, "QRFXY");
			this.DEEPRECHXY = GetNCVar(file, "DEEPRECHXY");
			this.QSPRINGXY = GetNCVar(file, "QSPRINGXY");
			this.QSLATXY = GetNCVar(file, "QSLATXY");
			this.QRFSXY = GetNCVar(file, "QRFSXY");
			this.QSPRINGSXY = GetNCVar(file, "QSPRINGSXY");
			this.RECHXY = GetNCVar(file, "RECHXY");
			this.FDEPTHXY = GetNCVar(file, "FDEPTHXY");
			this.RIVERCONDXY = GetNCVar(file, "RIVERCONDXY");
			this.RIVERBEDXY = GetNCVar(file, "RIVERBEDXY");
			this.EQZWT = GetNCVar(file, "EQZWT");
			this.PEXPXY = GetNCVar(file, "PEXPXY");
			this.ZLVLXY = GetNCVar(file, "ZLVLXY"); //2022年4月29日加，原Noah-MP没有此项。忽略该变量会导致restart后结果与正常执行有差距
			file.Close();
			
		}
		//		/// <summary>
		//		/// soil color: assuming a middle color category ?????????
		//		/// </summary>
		//		public static int ISC=4;
		//		/// <summary>
		//		///
		//		/// </summary>
		//		public static int IST=1;
		//		/// <summary>
		//		/// ice (sea ice = 1,land ice =-1, neither =0)
		//		/// </summary>
		//		public static int ICE;
		
		public void InitRestart()
		{
			double FIRA = 0;
			SMSTAV = new double[NRows, NCols];
			this.XICE = new double[NRows, NCols];
			SMSTOT = new double[NRows, NCols];
			if (TSLB == null)
				TSLB = new double[NRows, NSoil, NCols];
			for (int col = 0; col < NCols; col++) {

				int ITIMESTEP = 0;
				if (ITIMESTEP == 0) {
					for (int row = 0; row < NRows; row++) {
						if (SEAICE[row, col] > 0)
							XICE[row, col] = 1.0;
						
						if ((XLAND[row, col] - 1.5) >= 0) {    // Open water case
//							if (XICE[I, J] == 1 && IPRINT)
//								Console.WriteLine(" sea-ice at water point, I=" + I + "J=" + J);
							SMSTAV[row, col] = 1.0;
							SMSTOT[row, col] = 1.0;
							
							for (int K = 0; K < NSoil; K++) {
								SMOIS[row, K, col] = 1.0;
								TSLB[row, K, col] = 273.16;
							}
						} else {
							if (XICE[row, col] == 1.0) {        // Sea-ice case
								SMSTAV[row, col] = 1.0;
								SMSTOT[row, col] = 1.0;
								for (int K = 0; K < NSoil; K++) {
									SMOIS[row, K, col] = 1.0;
								}
							}
						}
					}
				}           // end of initialization over ocean


				//-----------------------------------------------------------------------

				//sea ice point ICE=1, land ice point ICE=-1, neither ICE=0;
				for (int row = 0; row < NRows; row++) {
					GridCell cell = Cells[row, col];
					//初始化土壤深度
					cell.ZSOIL = new FortDoubleArray(1, NSoil); //new double[NSoil];
					cell.ZSOIL[1] = -DZS[1];                    // depth to soil interfaces (<0) [m]
					for (int K = 2; K <= NSoil; K++) {
						
						cell.ZSOIL[K] = -DZS[K] + cell.ZSOIL[K - 1];
					}

					if (XICE[row, col] >= XICE_THRES) {
						cell.ICE = 1;                 // Sea-ice point
					} else if (IVGTYP[row, col] == ISICE)
						cell.ICE = -1;                // Land-ice point
					else
						cell.ICE = 0;                 // Neither sea ice or land ice.
					

//					if ((XLAND[row, col] - 1.5) >= 0)
//						continue;    	// Open water case

//					if (cell.ICE == 1) {
//
//						for (int K = 0; K < NSoil; K++) {
//							SH2O[row, K, col] = 1.0;
//						}
//						LAIXY[row, col] = 0.01;
//
//						continue; // Skip any processing at sea-ice points
//					} else {//     2D to 1D

					// IN only

					cell.COSZ = CALC_DECLIN(Driver.time0, latitude[row, col], longitude[row, col]);
					//cos zenith angle []
					
					cell.LATITUDE = latitude[row, col];                      // latitude [rad]
					cell.LONGITUDE = longitude[row, col];
					

					cell.ZLVL = ZLVL; //DZ8W[I, 1, J];                       // DZ8W: thickness of full levels; ZLVL forcing height [m]
					cell.VEGTYP = IVGTYP[row, col];                         // vegetation type
					cell.SOILTYP = ISLTYP[row, col];                          // soil type
					cell.FVEG = VEGFRA[row, col] / 100.0;                     // vegetation fraction [0-1]
					//cell.SHDMAX = SHDMAX[row, col] / 100.0;                   // Vegetation fraction annual max [0-1]
					cell.TBOT = TMN[row, col];                            // Fixed deep soil temperature for land
					
					//cell.SFCTMP = T3D[I, J];                       // temperature defined at intermediate level [K]
					
					//cell.Q2 = QV3D[I, J];
					;        // convert from mixing ratio to specific humidity [kg/kg]
					//cell.UU = U_PHY[I, J];                        // u-wind at interface [m/s]
					//cell.VV = V_PHY[I, J];                         // v-wind at interface [m/s]
					//cell.SOLDN = SWDOWN[I, J];                        // shortwave down from SW scheme [W/m2]
					//cell.LWDN = GLW[I, J];                           // total longwave down from LW scheme [W/m2]
					//cell.SFCPRS = P8W3D[I, J];  // surface pressure defined at intermediate level [Pa]
					//    consistent with temperature, mixing ratio
					//cell.PSFC = P8W3D[I, J];                        // surface pressure defined a full levels [Pa]
					//cell.PRCP = RAINBL[I, J] / DT;                      // timestep precipitation [mm/s]

					// IN/OUT fields
					//从nc中读入的NSnow范围为-2,-1,0，因此要加上Nsnow-1,转化为0,1,2
					//此处设为2，相当于原NoahMP中设为0。当原NoahMP中设为-2时，此处才能设为0
					cell.ISNOW = ISNOWXY[row, col];                // SNOW layers []
					
					cell.SH2O = new FortDoubleArray(1, NSoil);
					cell.STC = new FortDoubleArray(-NSnow + 1, NSoil);// new double[NSnow + NSoil];
					cell.SMC = new FortDoubleArray(1, NSoil);
					for (int z = 1; z <= NSoil; z++) {
						cell.SMC[z] = SMOIS[row, z - 1, col];  // soil total moisture [m3/m3]
						cell.SH2O[z] = SH2O[row, z - 1, col];  // soil liquid moisture [m3/m3]
						cell.STC[z] = TSLB[row, z - 1, col];  // soil temperatures [K]
						if (cell.SH2O[z] < 0) {
							throw new Exception("");
						}
						if (cell.STC[z] < 0 || cell.SH2O[z] < 0) {
							cell.NeedCal = false;
							cell.TV = 273;
							cell.STC[z] = 273;
							cell.SH2O[z] = 0.23;
							cell.SMC[z] = 0.23;
						}
					}
					for (int z = 1 - NSnow; z <= 0; z++) {
						if (TSNOXY != null)
							cell.STC[z] = TSNOXY[row, z + (NSnow - 1), col]; // SNOW temperatures [K]
						else
							cell.STC[z] = 273;
						if (cell.STC[z] < 100)
							cell.STC[z] = 273;
						//throw new Exception("");
					}
					cell.SNEQV = SNOW[row, col];             // SNOW water equivalent [mm]
					cell.SNOWH = SNOWH[row, col];               // SNOW depth [m]
					cell.QSFC = QSFC[row, col];
					
					
					//Console.WriteLine("QSFC:"+cell.QSFC);
					if (double.IsNaN(cell.QSFC) || double.IsInfinity(cell.QSFC)) {
						cell.QSFC = 0.005;
					}
					// INOUT (with no Noah LSM equivalent)
					if (TVXY != null)
						cell.TV = TVXY[row, col];                // leaf temperature [K]
					else
						cell.TV = 273;
					if (TGXY != null) {
						cell.TG = TGXY[row, col];               // ground temperature [K]
						if (double.IsNaN(cell.TG))
							cell.TG = 273 + 10;
					} else {
						cell.TG = 273;
					}
//						if(double. cell.TG)
//							throw new Exception("");
					cell.CANLIQ = CANLIQXY[row, col];              // canopy liquid water [mm]
					cell.CANICE = CANICEXY[row, col];              // canopy frozen water [mm]
					cell.EAH = EAHXY[row, col];               // canopy vapor pressure [Pa]
					cell.TAH = TAHXY[row, col];              // canopy temperature [K]
					cell.CM = CMXY[row, col];              // avg. momentum exchange (MP only) [m/s]
					cell.CH = CHXY[row, col];              // avg. heat exchange (MP only) [m/s]
					cell.FWET = FWETXY[row, col];             // canopy fraction wet or SNOW
					cell.SNEQVO = SNEQVOXY[row, col];             // SWE previous timestep
					cell.ALBOLD = ALBOLDXY[row, col];               // albedo previous timestep, for SNOW aging
					cell.QSNOW = QSNOWXY[row, col];             // SNOW falling on ground
					cell.WSLAKE = WSLAKEXY[row, col];              // lake water storage (can be neg.) (mm)
					cell.ZWT = ZWTXY[row, col];               // depth to water table [m]
					cell.WA = WAXY[row, col];               // water storage in aquifer [mm]
					cell.WT = WTXY[row, col];               // water in aquifersaturated soil [mm]

					double sum = 0;
					cell.ZSNSO = new FortDoubleArray(1 - NSnow, NSoil); //new double[NSnow + NSoil];
					for (int z = 1 - NSnow; z <= NSoil; z++) {
						cell.ZSNSO[z] = ZSNSOXY[row, z + (NSnow - 1), col]; // depth to layer interface
						if (cell.ZSNSO[z] < -1000) {
							cell.ZSNSO.data = new double[]{ 0, 0, 0, 0.1, 0.4, 1.0, 2 };
							break;
						}
					}
					
					
					cell.SNICE = new FortDoubleArray(1 - NSnow, 0);// new double[NSnow];
					//Console.WriteLine("SNICE printed "+cell.ICE);
					cell.SNLIQ = new FortDoubleArray(1 - NSnow, 0);// new double[NSnow];
					for (int z = 1 - NSnow; z <= 0; z++) {
						cell.SNICE[z] = SNICEXY[row, z + (NSnow - 1), col];// SNOW layer ice content
						cell.SNLIQ[z] = SNLIQXY[row, z + (NSnow - 1), col]; // SNOW layer water content
						if (cell.SNICE[z] < -1e10) {
							cell.SNICE[z] = 0;
							cell.SNLIQ[z] = 0;
						}
					}
					cell.LFMASS = LFMASSXY[row, col];               // leaf mass
					cell.RTMASS = RTMASSXY[row, col];               // root mass
					cell.STMASS = STMASSXY[row, col];              // stem mass
					cell.WOOD = WOODXY[row, col];               // mass of wood (incl. woody roots) [g/m2]
					cell.STBLCP = STBLCPXY[row, col];              // stable carbon pool
					cell.FASTCP = FASTCPXY[row, col];                // fast carbon pool
					cell.LAI = LAIXY[row, col];                // leaf area index [-] (no SNOW effects)
					cell.SAI = XSAIXY[row, col];                // stem area index [-] (no SNOW effects)
					cell.TAUSS = TAUSSXY[row, col];                // non-dimensional SNOW age
					if (NoahMP.OPT_RUN == 5) {
						cell.SMCEQ = new FortDoubleArray(1, NSoil);// new double[NSoil];
						for (int z = 0; z < NSoil; z++) {
							cell.SMCEQ[z] = SMOISEQ[row, z, col];
						}
						cell.SMCWTD = SMCWTDXY[row, col];
					}
//						double RECH = 0;
//						double DEEPRECH = 0;

					// Initialized local

					cell.FICEOLD = new FortDoubleArray(1 - NSnow, 0);// new double[NSnow];
					for (int z = cell.ISNOW + 1; z <= 0; z++) {
						cell.FICEOLD[z] = SNICEXY[row, z + (NSnow - 1), col] / (SNICEXY[row, z + (NSnow - 1), col] + SNLIQXY[row, z + (NSnow - 1), col]);   // SNOW ice fraction
//						if(double.IsNaN(cell.FICEOLD[z]))
//							throw new Exception("");
						
					}
					cell.IST = 1;                                   // MP surface type: 1 = land; 2 = lake
					cell.ISC = 4;                                   // soil color: assuming a middle color category ?????????


					if (cell.SOILTYP == 14 && XICE[row, col] == 0) {
//							if (IPRINT)
//								Console.WriteLine(" SOIL TYPE FOUND TO BE WATER AT A LAND-POINT");
//							if (IPRINT)
//								Console.WriteLine(I + " " + J + "RESET SOIL in surfce.F");
						cell.SOILTYP = 7;
					}

					if (IVGTYP[row, col] == ISURBAN || IVGTYP[row, col] == 31 ||
					    IVGTYP[row, col] == 32 || IVGTYP[row, col] == 33) {
						cell.VEGTYP = ISURBAN;
					}
					
					if (cell.VEGTYP == 25)
						cell.FVEG = 0.0;                  // Set playa, lava, sand to bare
					if (cell.VEGTYP == 25)
						cell.LAI = 0.0;
					if (cell.VEGTYP == 26)
						cell.FVEG = 0.0;                 // hard coded for USGS
					if (cell.VEGTYP == 26)
						cell.LAI = 0.0;
					if (cell.VEGTYP == 27)
						cell.FVEG = 0.0;
					if (cell.VEGTYP == 27)
						cell.LAI = 0.0;

				} // } of land-sea test

			}// ENDDO ILOOP                                                       // of I loop
			//}//  ENDDO JLOOP
//			StreamWriter sw = new StreamWriter("d:/liu.txt");
//			for (int row = 0; row < NRows; row++) {
//				for (int col = 0; col < NCols; col++) {
//					if (Cells[row, col].ZSNSO != null)
//						sw.Write(Cells[row, col].SH2O[1] + " ");
//					else
//						sw.Write("-9999 ");
//				}
//				sw.WriteLine();
//			}
//			sw.Close();
			
		}
		public void AssignToCell_Restart()
		{
			for (int col = 0; col < NCols; col++) {
				for (int row = 0; row < NRows; row++) {
					GridCell cell = Cells[row, col];
					//处理土壤温度
					cell.STC = new FortDoubleArray(1 - NSnow, NSoil);
					for (int z = 1; z <= NSoil; z++) {
						cell.STC[z] = TSLB[row, z - 1, col];
					}
					//this.TSLB = GetNCVar3D_Restart(file, "SOIL_T");
					//处理雪温
					for (int z = 1 - NSnow; z <= 0; z++) {
						cell.STC[z] = TSNOXY[row, z + (NSnow - 1), col];
					}
					//this.TSNOXY = GetNCVar3D_Restart(file, "SNOW_T");
					//处理土壤湿度
					cell.SMC = new FortDoubleArray(1, NSoil);
					for (int z = 1; z <= NSoil; z++) {
						cell.SMC[z] = SMOIS[row, z - 1, col];
					}
					//this.SMOIS = GetNCVar3D_Restart(file, "SMC");
					//液态土壤湿度，猜测可能是不包括冰在内
					cell.SH2O = new FortDoubleArray(1, NSoil);
					for (int z = 1; z <= NSoil; z++) {
						cell.SH2O[z] = SH2O[row, z - 1, col];
						if (cell.SH2O[z] < 0)
							throw new Exception("看看SH2O的数据是不是有问题");
					}
					//this.SH2O = GetNCVar3D_Restart(file, "SH2O");
					//雪和土壤层的深度
					cell.ZSNSO = new FortDoubleArray(1 - NSnow, NSoil);
					for (int z = 1 - NSnow; z <= NSoil; z++) {
						cell.ZSNSO[z] = ZSNSOXY[row, z + (NSnow - 1), col];
					}
					//this.ZSNSOXY = GetNCVar3D_Restart(file, "ZSNSO");
					//雪层中冰厚度
					for (int z = 1 - NSnow; z <= 0; z++) {
						cell.SNICE[z] = SNICEXY[row, z + (NSnow - 1), col];
					}
					//this.SNICEXY = GetNCVar3D_Restart(file, "SNICE");
					//雪中液态厚度
					for (int z = 1 - NSnow; z <= 0; z++) {
						cell.SNLIQ[z] = SNLIQXY[row, z + (NSnow - 1), col];
					}
					cell.HGT = TERRAIN[row, col];
					//this.SNLIQXY = GetNCVar3D_Restart(file, "SNLIQ");
					cell.QSNOW = QSNOWXY[row, col];    //this.QSNOWXY = GetNCVar(file, "QSNOW");
					cell.FWET = FWETXY[row, col];   	//this.FWETXY = GetNCVar(file, "FWET");
					cell.SNEQVO = SNEQVOXY[row, col]; //this.SNEQVOXY = GetNCVar(file, "SNEQVO");
					cell.EAH = EAHXY[row, col];        //this.EAHXY = GetNCVar(file, "EAH");
					cell.TAH = TAHXY[row, col]; //this.TAHXY = GetNCVar(file, "TAH");
					cell.ALBOLD = ALBOLDXY[row, col];	//this.ALBOLDXY = GetNCVar(file, "ALBOLD");
					cell.CM = CMXY[row, col];	//this.CMXY = GetNCVar(file, "CM");
					cell.CH = CHXY[row, col]; //this.CHXY = GetNCVar(file, "CH");
					cell.ISNOW = ISNOWXY[row, col]; //this.ISNOWXY = GetNCVar_Int(file, "ISNOW");
					cell.CANLIQ = CANLIQXY[row, col];	//this.CANLIQXY = GetNCVar(file, "CANLIQ");
					cell.CANICE = CANICEXY[row, col];	//this.CANICEXY = GetNCVar(file, "CANICE");
					cell.SNEQV = SNOW[row, col]; //this.SNOW = GetNCVar(file, "SNEQV");
					cell.SNOWH = SNOWH[row, col];//this.SNOWH = GetNCVar(file, "SNOWH");
					cell.TV = TVXY[row, col];	//this.TVXY = GetNCVar(file, "TV");
					cell.TG = TGXY[row, col]; //this.TGXY = GetNCVar(file, "TG");
					cell.ZWT = ZWTXY[row, col]; //this.ZWTXY = GetNCVar(file, "ZWT");
					cell.WA = WAXY[row, col]; //this.WAXY = GetNCVar(file, "WA");
					cell.WT = WTXY[row, col];//this.WTXY = GetNCVar(file, "WT");
					cell.WSLAKE = WSLAKEXY[row, col];//this.WSLAKEXY = GetNCVar(file, "WSLAKE");
					cell.LFMASS = LFMASSXY[row, col];//this.LFMASSXY = GetNCVar(file, "LFMASS");
					cell.RTMASS = RTMASSXY[row, col];//this.RTMASSXY = GetNCVar(file, "RTMASS");
					cell.STMASS = STMASSXY[row, col];//this.STMASSXY = GetNCVar(file, "STMASS");
					cell.WOOD = WOODXY[row, col];//this.WOODXY = GetNCVar(file, "WOOD");
					cell.STBLCP = STBLCPXY[row, col];//this.STBLCPXY = GetNCVar(file, "STBLCP");
					cell.FASTCP = FASTCPXY[row, col]; //this.FASTCPXY = GetNCVar(file, "FASTCP");
					cell.LAI = LAIXY[row, col]; //this.LAIXY = GetNCVar(file, "LAI");
					cell.SAI = XSAIXY[row, col]; //this.XSAIXY = GetNCVar(file, "SAI");
					cell.FVEG = VEGFRA[row, col] * 0.01; //this.VEGFRA = GetNCVar(file, "FPAR");
					cell.FOLN = 1.0;
					//GVFMIN好像并不参与计算，暂可以不考虑 //this.GVFMIN = GetNCVar(file, "GVFMIN");
					cell.SHDMAX = SHDMAX[row, col] / 100;//this.SHDMAX = GetNCVar(file, "SHDMAX");
					if (SHDMIN != null) {
						cell.SHDMIN = SHDMIN[row, col] / 100; //this.SHDMIN = GetNCVar(file, "SHDMIN");
					}
					//GVFMAX好像也不参与计算，暂可以不考虑; //this.GVFMAX = GetNCVar(file, "GVFMAX");
					//累积量应该不用更新进去,ACSNOM[row,col]; //this.ACSNOM = GetNCVar(file, "ACMELT");
					//this.ACSNOW = GetNCVar(file, "ACSNOW");
					cell.TAUSS = TAUSSXY[row, col];//this.TAUSSXY = GetNCVar(file, "TAUSS");
					cell.QSFC = QSFC[row, col]; //this.QSFC = GetNCVar(file, "QSFC");
					//this.ACCSFCRUNOFF = GetNCVar(file, "SFCRUNOFF");
					//this.ACCUGDRUNOFF = GetNCVar(file, "UGDRUNOFF");
					//this.ACCPRCP = GetNCVar(file, "ACCPRCP");
					//this.ACCECAN = GetNCVar(file, "ACCECAN");
					//this.ACCEDIR = GetNCVar(file, "ACCEDIR");
					//this.ACCETRAN = GetNCVar(file, "ACCETRAN");
					
					//应该无实际意义  //this.SMOISEQ = GetNCVar3D_Restart(file, "SMOISEQ");
					//this.AREAXY = GetNCVar(file, "AREAXY");
					//this.SMCWTDXY = GetNCVar(file, "SMCWTDXY");
					//this.QRFXY = GetNCVar(file, "QRFXY");
					//this.DEEPRECHXY = GetNCVar(file, "DEEPRECHXY");
					cell.ZLVL = ZLVLXY[row, col]; //2022年4月29日加，原Noah-MP没有此项。忽略该变量会导致restart后结果与正常执行有差距
					if (NoahMP.OPT_RUN == 5) {
						//this.QSPRINGXY = GetNCVar(file, "QSPRINGXY");
						//this.QSLATXY = GetNCVar(file, "QSLATXY");
						//this.QRFSXY = GetNCVar(file, "QRFSXY");
						//this.QSPRINGSXY = GetNCVar(file, "QSPRINGSXY");
						//this.RECHXY = GetNCVar(file, "RECHXY");
						//this.FDEPTHXY = GetNCVar(file, "FDEPTHXY");
						//this.RIVERCONDXY = GetNCVar(file, "RIVERCONDXY");
						//this.RIVERBEDXY = GetNCVar(file, "RIVERBEDXY");
						//this.EQZWT = GetNCVar(file, "EQZWT");
						//this.PEXPXY = GetNCVar(file, "PEXPXY");
					}
				}
			}
		}
		public void AssignToCell_NoRestart()
		{
			double SLOPETYP = 1;                               // set underground runoff slope term
			double FIRA = 0;
			SMSTAV = new double[NRows, NCols];
			this.XICE = new double[NRows, NCols];
			SMSTOT = new double[NRows, NCols];
			if (TSLB == null)
				TSLB = new double[NRows, NSoil, NCols];
			for (int col = 0; col < NCols; col++) {
				int ITIMESTEP = 0;
				if (ITIMESTEP == 0) {
					for (int row = 0; row < NRows; row++) {
						if (SEAICE[row, col] > 0)
							XICE[row, col] = 1.0;
						if ((XLAND[row, col] - 1.5) >= 0) {    // Open water case
//							if (XICE[I, J] == 1 && IPRINT)
//								Console.WriteLine(" sea-ice at water point, I=" + I + "J=" + J);
							SMSTAV[row, col] = 1.0;
							SMSTOT[row, col] = 1.0;
							
							for (int K = 0; K < NSoil; K++) {
								SMOIS[row, K, col] = 1.0;
								TSLB[row, K, col] = 273.16;
							}
						} else {
							if (XICE[row, col] == 1.0) {        // Sea-ice case
								SMSTAV[row, col] = 1.0;
								SMSTOT[row, col] = 1.0;
								for (int K = 0; K < NSoil; K++) {
									SMOIS[row, K, col] = 1.0;
								}
							}
						}
					}
				}           // end of initialization over ocean


				//-----------------------------------------------------------------------

				//sea ice point ICE=1, land ice point ICE=-1, neither ICE=0;
				for (int row = 0; row < NRows; row++) {
					GridCell cell = Cells[row, col];
					//初始化土壤深度
					cell.ZSOIL = new FortDoubleArray(1, NSoil); //new double[NSoil];
					cell.ZSOIL[1] = -DZS[1];                    // depth to soil interfaces (<0) [m]
					for (int K = 2; K <= NSoil; K++) {
						cell.ZSOIL[K] = -DZS[K] + cell.ZSOIL[K - 1];
					}

					if (XICE[row, col] >= XICE_THRES) {
						cell.ICE = 1;                 // Sea-ice point
					} else if (IVGTYP[row, col] == ISICE)
						cell.ICE = -1;                // Land-ice point
					else
						cell.ICE = 0;                 // Neither sea ice or land ice.
					

//					if ((XLAND[row, col] - 1.5) >= 0)
//						continue;    	// Open water case

//					if (cell.ICE == 1) {
//
//						for (int K = 0; K < NSoil; K++) {
//							SH2O[row, K, col] = 1.0;
//						}
//						LAIXY[row, col] = 0.01;
//
//						continue; // Skip any processing at sea-ice points
//					} else {//     2D to 1D

					// IN only

					cell.COSZ = CALC_DECLIN(Driver.time0, latitude[row, col], longitude[row, col]);
					//cos zenith angle []
					
					cell.LATITUDE = latitude[row, col];                      // latitude [rad]
					cell.LONGITUDE = longitude[row, col];
					cell.HGT = TERRAIN[row, col];
					cell.ZLVL = ZLVL; //DZ8W[I, 1, J];                       // DZ8W: thickness of full levels; ZLVL forcing height [m]
					cell.VEGTYP = IVGTYP[row, col];                         // vegetation type
					cell.SOILTYP = ISLTYP[row, col];                          // soil type
					cell.FVEG = VEGFRA[row, col] / 100.0;                     // vegetation fraction [0-1]
					cell.SHDMAX = SHDMAX[row, col] / 100;            // Vegetation fraction annual max [0-1]
					cell.SHDMIN = SHDMIN[row, col] / 100;
					cell.TBOT = TMN[row, col];                            // Fixed deep soil temperature for land
					
					//cell.SFCTMP = T3D[I, J];                       // temperature defined at intermediate level [K]
					
					//cell.Q2 = QV3D[I, J];
					;        // convert from mixing ratio to specific humidity [kg/kg]
					//cell.UU = U_PHY[I, J];                        // u-wind at interface [m/s]
					//cell.VV = V_PHY[I, J];                         // v-wind at interface [m/s]
					//cell.SOLDN = SWDOWN[I, J];                        // shortwave down from SW scheme [W/m2]
					//cell.LWDN = GLW[I, J];                           // total longwave down from LW scheme [W/m2]
					//cell.SFCPRS = P8W3D[I, J];  // surface pressure defined at intermediate level [Pa]
					//    consistent with temperature, mixing ratio
					//cell.PSFC = P8W3D[I, J];                        // surface pressure defined a full levels [Pa]
					//cell.PRCP = RAINBL[I, J] / DT;                      // timestep precipitation [mm/s]

					// IN/OUT fields
					//从nc中读入的NSnow范围为-2,-1,0，因此要加上Nsnow,转化为0,1,2
					cell.ISNOW = ISNOWXY[row, col];                // SNOW layers []					
					
					cell.STC = new FortDoubleArray(1 - NSnow, NSoil);  //第3,4,5,6是土壤层
					cell.SMC = new FortDoubleArray(1, 4); //new double[NSoil];
					cell.SH2O = new FortDoubleArray(1, 4); //new double[NSoil];
					for (int z = 1; z <= NSoil; z++) {
						cell.SMC[z] = SMOIS[row, z - 1, col];  // soil total moisture [m3/m3]
						cell.SH2O[z] = SH2O[row, z - 1, col];  // soil liquid moisture [m3/m3]
						cell.STC[z] = TSLB[row, z - 1, col];  // soil temperatures [K]
						if (cell.SH2O[z] < 1e-10) {
							cell.SH2O[z] = cell.SMC[z];
						}
						if (Math.Abs(cell.SH2O[z]) > 500) {
							cell.SH2O[z] = 0.25;
						}
						if (cell.SH2O[z] < 0)
							throw new Exception();
					}
//					for (int z = 0; z < NSnow; z++) {
//
//						// SNOW temperatures [K]
//
//						cell.STC[z] = TSLB[row,z,col];
//					}
					cell.SNEQV = SNOW[row, col];             // SNOW water equivalent [mm]
					cell.SNOWH = SNOWH[row, col];               // SNOW depth [m]
					cell.QSFC = QSFC[row, col];
					cell.PSFC = PSFC[row, col];
					//Console.WriteLine("QSFC:"+cell.QSFC);
					if (double.IsNaN(cell.QSFC) || double.IsInfinity(cell.QSFC)) {
						cell.QSFC = 0.005;
					}
					// INOUT (with no Noah LSM equivalent)
					if (TVXY != null)
						cell.TV = TVXY[row, col];                // leaf temperature [K]
					else
						cell.TV = 273;
					
					cell.TG = TGXY[row, col];					
					cell.CANLIQ = CANLIQXY[row, col];              // canopy liquid water [mm]
					cell.CANICE = CANICEXY[row, col];              // canopy frozen water [mm]
					cell.EAH = EAHXY[row, col];               // canopy vapor pressure [Pa]
					cell.TAH = TAHXY[row, col];              // canopy temperature [K]
					cell.CM = CMXY[row, col];              // avg. momentum exchange (MP only) [m/s]
					cell.CH = CHXY[row, col];              // avg. heat exchange (MP only) [m/s]
					cell.FWET = FWETXY[row, col];             // canopy fraction wet or SNOW
					cell.SNEQVO = SNEQVOXY[row, col];             // SWE previous timestep
					cell.ALBOLD = ALBOLDXY[row, col];               // albedo previous timestep, for SNOW aging
					cell.QSNOW = QSNOWXY[row, col];             // SNOW falling on ground
					cell.WSLAKE = WSLAKEXY[row, col];              // lake water storage (can be neg.) (mm)
					cell.ZWT = ZWTXY[row, col];               // depth to water table [m]
					cell.WA = WAXY[row, col];               // water storage in aquifer [mm]
					cell.WT = WTXY[row, col];               // water in aquifersaturated soil [mm]

					double sum = 0;
					cell.ZSNSO = new FortDoubleArray(1 - NSnow, NSoil);// new double[NSnow + NSoil];
					cell.ZSNSO.data = new double[]{ 0, 0, 0, -0.1, -0.4, -1.0, -2 };
					
					
					
					cell.SNICE = new FortDoubleArray(1 - Driver.NSnow, 0);// new double[NSnow];
					//Console.WriteLine("SNICE printed "+cell.ICE);
					cell.SNLIQ = new FortDoubleArray(1 - Driver.NSnow, 0);//= new double[NSnow];
//					for (int z = 0; z < NSnow; z++) {
//						cell.SNICE[z] = SNICEXY[row, z, col];// SNOW layer ice content
//						cell.SNLIQ[z] = SNLIQXY[row, z, col]; // SNOW layer water content
//						if (cell.SNICE[z] < -1e10) {
//							cell.SNICE[z] = 0;
//							cell.SNLIQ[z] = 0;
//						}
//					}
					//cell.LFMASS = LFMASSXY[row, col];               // leaf mass
					//cell.RTMASS = RTMASSXY[row, col];               // root mass
					//cell.STMASS = STMASSXY[row, col];              // stem mass
					//cell.WOOD = WOODXY[row, col];               // mass of wood (incl. woody roots) [g/m2]
					//cell.STBLCP = STBLCPXY[row, col];              // stable carbon pool
					//cell.FASTCP = FASTCPXY[row, col];                // fast carbon pool
					cell.LAI = LAIXY[row, col];                // leaf area index [-] (no SNOW effects)
					//cell.SAI = XSAIXY[row, col];                // stem area index [-] (no SNOW effects)
					//cell.TAUSS = TAUSSXY[row, col];                // non-dimensional SNOW age
					cell.SMCEQ = new FortDoubleArray(1, NSoil);
					for (int z = 1; z <= NSoil; z++) {
						cell.SMCEQ[z] = 0; //SMOISEQ[row, z, col];
					}
//					cell.SMCWTD = SMCWTDXY[row, col];
//						double RECH = 0;
//						double DEEPRECH = 0;

					// Initialized local

					cell.FICEOLD = new FortDoubleArray(1 - NSnow, 0);// new double[NSnow];
					
					double CO2 = 395e-06;
					double O2 = 0.209;
					double CO2PP = CO2 * cell.SFCPRS;                         // partial pressure co2 [Pa]
					double O2PP = O2 * cell.SFCPRS;                          // partial pressure  o2 [Pa]
					cell.FOLN = 1.0;                              // for now, set to nitrogen saturation
					double QC = -1e36;                     // test dummy value
					double PBLH = 0;                     // test dummy value // PBL height
					double DZ8W1D = 2 * cell.ZLVL;                          // thickness of atmospheric layers
					cell.SLOPETYP = 2;                               // set underground runoff slope term
					cell.IST = 1;                                   // MP surface type: 1 = land; 2 = lake
					cell.ISC = 4;                                   // soil color: assuming a middle color category ?????????


					if (cell.SOILTYP == 14 && XICE[row, col] == 0) {
//							if (IPRINT)
//								Console.WriteLine(" SOIL TYPE FOUND TO BE WATER AT A LAND-POINT");
//							if (IPRINT)
//								Console.WriteLine(I + " " + J + "RESET SOIL in surfce.F");
						cell.SOILTYP = 7;
					}

					if (IVGTYP[row, col] == ISURBAN || IVGTYP[row, col] == 31 ||
					    IVGTYP[row, col] == 32 || IVGTYP[row, col] == 33) {
						cell.VEGTYP = vegparams.ISURBAN;
					}
					
					if (cell.VEGTYP == 25)
						cell.FVEG = 0.0;                  // Set playa, lava, sand to bare
					if (cell.VEGTYP == 25)
						cell.LAI = 0.0;
					if (cell.VEGTYP == 26)
						cell.FVEG = 0.0;                 // hard coded for USGS
					if (cell.VEGTYP == 26)
						cell.LAI = 0.0;
					if (cell.VEGTYP == 27)
						cell.FVEG = 0.0;
					if (cell.VEGTYP == 27)
						cell.LAI = 0.0;

				} // } of land-sea test

			}// ENDDO ILOOP                                                       // of I loop
			//}//  ENDDO JLOOP
			
			
		}
		int outputMode = 0;
		public void add_to_output(int ncid, double[,] data, string name, string description, string unit)
		{
			if (!outputParams.Contains(name))
				return;
			if (data == null) {
				Console.WriteLine("这个二维变量不存在" + name);
				return;
			}
			if (outputMode == 0) {
				int varid = 0;
				NC_LowAPI.nc_def_var(ncid, name, NC_Type.NC_FLOAT, 2, new int[]{ 0, 1 }, ref varid);
				NC_LowAPI.nc_put_att_string(ncid, varid, "description", description.Length, description);
				NC_LowAPI.nc_put_att_string(ncid, varid, "unit", unit.Length, unit);
			} else {
				int varid = -1;
				NC_LowAPI.nc_inq_varid(ncid, name, ref varid);
				float[] array = new float[NRows * NCols];
				for (int row = 0; row < NRows; row++) {
					for (int col = 0; col < NCols; col++) {
						array[row * NCols + col] = (float)data[row, col];
					}
				}
				int error = NC_LowAPI.nc_put_var_float(ncid, varid, array);
			}
			
			
			//需要补充属性的赋值
			
			//NC_LowAPI.nc_close(ncid);
		}
		public void add_to_output(int ncid, int[,] data, string name, string description, string unit)
		{
			if (!outputParams.Contains(name))
				return;
			if (data == null) {
				Console.WriteLine("这个二维变量不存在" + name);
				return;
			}
			if (outputMode == 0) {
				int varid = 0;
				NC_LowAPI.nc_def_var(ncid, name, NC_Type.NC_INT, 2, new int[]{ 0, 1 }, ref varid);
				NC_LowAPI.nc_put_att_string(ncid, varid, "description", description.Length, description);
				NC_LowAPI.nc_put_att_string(ncid, varid, "unit", unit.Length, unit);
			} else {
				int varid = 0;
				NC_LowAPI.nc_inq_varid(ncid, name, ref varid);
				int[] array = new int[NRows * NCols];
				for (int row = 0; row < NRows; row++) {
					for (int col = 0; col < NCols; col++) {
						array[row * NCols + col] = data[row, col];
					}
				}
				NC_LowAPI.nc_put_var_int(ncid, varid, array);
			}
			
			//需要补充属性的赋值
			
			//NC_LowAPI.nc_close(ncid);
		}
		public void add_to_output_3d(int ncid, int[,,] data, int[] dimIDs, int NLayers, string name, string description, string unit)
		{
			if (!outputParams.Contains(name))
				return;
			if (data == null) {
				Console.WriteLine("这个二维变量不存在" + name);
				return;
			}
			if (outputMode == 0) {
				int varid = 0;
				NC_LowAPI.nc_def_var(ncid, name, NC_Type.NC_INT, 3, dimIDs, ref varid);
				NC_LowAPI.nc_put_att_string(ncid, varid, "description", description.Length, description);
				NC_LowAPI.nc_put_att_string(ncid, varid, "unit", unit.Length, unit);
			} else {
				int varid = 0;
				NC_LowAPI.nc_inq_varid(ncid, name, ref varid);
				int[] array = new int[data.Length];
				for (int layer = 0; layer < NLayers; layer++) {
					for (int row = 0; row < NRows; row++) {
						for (int col = 0; col < NCols; col++) {
							array[layer * NRows * NCols + row * NCols + col] = data[row, layer, col];
						}
					}
				}
				NC_LowAPI.nc_put_var_int(ncid, varid, array);
			}
			//需要补充属性的赋值
		}
		public void add_to_output_3d(int ncid, double[,,] data, int[] dimIDs, int NLayers, string name, string description, string unit)
		{
			if (!outputParams.Contains(name))
				return;
			if (data == null) {
				Console.WriteLine("这个二维变量不存在" + name);
				return;
			}
			if (outputMode == 0) {
				int varid = 0;
				NC_LowAPI.nc_def_var(ncid, name, NC_Type.NC_FLOAT, 3, dimIDs, ref varid);
				NC_LowAPI.nc_put_att_string(ncid, varid, "description", description.Length, description);
				NC_LowAPI.nc_put_att_string(ncid, varid, "unit", unit.Length, unit);
			} else {
				int varid = 0;
				NC_LowAPI.nc_inq_varid(ncid, name, ref varid);
				float[] array = new float[data.Length];
				for (int layer = 0; layer < NLayers; layer++) {
					for (int row = 0; row < NRows; row++) {
						for (int col = 0; col < NCols; col++) {
							//array[row * NLayers * NCols + layer * NCols + col] = (float)data[row, layer, col];
							array[layer * NRows * NCols + row * NCols + col] = (float)data[row, layer, col];
						}
					}
				}
				NC_LowAPI.nc_put_var_float(ncid, varid, array);
			}
			//需要补充属性的赋值
		}
		public void add_to_restart(int ncid, int[,] data, int[] dimIDs, string name, string description, string unit)
		{
			if (outputMode == 0) {
				int varid = 0;
				NC_LowAPI.nc_def_var(ncid, name, NC_Type.NC_INT, 2, dimIDs, ref varid);
				NC_LowAPI.nc_put_att_string(ncid, varid, "description", description.Length, description);
				NC_LowAPI.nc_put_att_string(ncid, varid, "unit", unit.Length, unit);
			} else {
				int varid = 0;
				NC_LowAPI.nc_inq_varid(ncid, name, ref varid);
				int[] array = new int[NRows * NCols];
				for (int row = 0; row < NRows; row++) {
					for (int col = 0; col < NCols; col++) {
						array[row * NCols + col] = data[row, col];
					}
				}
				NC_LowAPI.nc_put_var_int(ncid, varid, array);
			}
		}
		public void add_to_restart(int ncid, double[,] data, int[] dimIDs, string name, string description, string unit)
		{
			if (data == null) {
				return;
			}
			if (outputMode == 0) {
				int varid = 0;
				NC_LowAPI.nc_def_var(ncid, name, NC_Type.NC_FLOAT, 2, dimIDs, ref varid);
				NC_LowAPI.nc_put_att_string(ncid, varid, "description", description.Length, description);
				NC_LowAPI.nc_put_att_string(ncid, varid, "unit", unit.Length, unit);
			} else {
				int varid = 0;
				NC_LowAPI.nc_inq_varid(ncid, name, ref varid);
				float[] array = new float[NRows * NCols];
				for (int row = 0; row < NRows; row++) {
					for (int col = 0; col < NCols; col++) {
						array[row * NCols + col] = (float)data[row, col];
					}
				}
				NC_LowAPI.nc_put_var_float(ncid, varid, array);
			}
		}
		public void add_to_restart_3d(int ncid, int[,,] data, int[] dimIDs, int NLayers, string name, string description, string unit)
		{
			if (data == null) {
				return;
			}
			if (outputMode == 0) {
				int varid = 0;
				NC_LowAPI.nc_def_var(ncid, name, NC_Type.NC_INT, 3, dimIDs, ref varid);
				NC_LowAPI.nc_put_att_string(ncid, varid, "description", description.Length, description);
				NC_LowAPI.nc_put_att_string(ncid, varid, "unit", unit.Length, unit);
			} else {
				int varid = 0;
				NC_LowAPI.nc_inq_varid(ncid, name, ref varid);
				
				float[] array = new float[data.Length];
				for (int layer = 0; layer < NLayers; layer++) {
					for (int row = 0; row < NRows; row++) {
						for (int col = 0; col < NCols; col++) {
							array[row * NLayers * NCols + layer * NCols + col] = (float)data[row, layer, col];
						}
					}
				}
				NC_LowAPI.nc_put_var_float(ncid, varid, array);
			}
		}
		public void add_to_restart_3d(int ncid, double[,,] data, int[] dimIDs, int NLayers, string name, string description, string unit)
		{
			if (data == null) {
				return;
			}
			if (outputMode == 0) {
				int varid = 0;
				NC_LowAPI.nc_def_var(ncid, name, NC_Type.NC_FLOAT, 3, dimIDs, ref varid);
				NC_LowAPI.nc_put_att_string(ncid, varid, "description", description.Length, description);
				NC_LowAPI.nc_put_att_string(ncid, varid, "unit", unit.Length, unit);
			} else {
				int varid = -1;
				NC_LowAPI.nc_inq_varid(ncid, name, ref varid);
				float[] array = new float[data.Length];
				for (int layer = 0; layer < NLayers; layer++) {
					for (int row = 0; row < NRows; row++) {
						for (int col = 0; col < NCols; col++) {
							array[row * NLayers * NCols + layer * NCols + col] = (float)data[row, layer, col];
						}
					}
				}
				NC_LowAPI.nc_put_var_float(ncid, varid, array);
			}
			
			//需要补充属性的赋值
			
			//NC_LowAPI.nc_close(ncid);
		}
		double[,] TSK;
		double[,] HFX;
		double[,] GRDFLX;
		double[,] ALBEDO;
		double[,] SNOWC;
		double[,] CANWAT;
		double[,] EMISS;
		double[,] FSAXY;
		double[,] FIRAXY;
		double[,] ECANXY;
		double[,]	ETRANXY;
		double[,] EDIRXY;
		
		double[,] T2MVXY;
		double[,] T2MBXY;
		double[,]	Q2MVXY;
		// specific humidity to mixing ratio
		double[,] Q2MBXY;
		// consistent with registry def of Q2
		double[,] TRADXY;
		double[,] NEEXY;
		double[,] GPPXY;
		double[,] NPPXY;
		double[,] FVEGXY;
		double[,]	RUNSFXY;
		double[,] RUNSBXY;
		//double[,]			RUNSBXY =new double[NRows,NCols];
		double[,]	APARXY;
		double[,] PSNXY;
		double[,] SAVXY;
		double[,] SAGXY;
		double[,] RSSUNXY;
		double[,] RSSHAXY;
		double[,] BGAPXY;
		double[,] WGAPXY;
		double[,] TGVXY;
		double[,] TGBXY;
		double[,] CHVXY;
		double[,] CHBXY;
		double[,] IRCXY;
		double[,] IRGXY;
		double[,] SHCXY;
		double[,] SHGXY;
		double[,] EVGXY;
		double[,] GHVXY;
		double[,] IRBXY;
		double[,] SHBXY;
		double[,] EVBXY;
		double[,] GHBXY;
		double[,] TRXY;
		double[,] EVCXY;
		double[,] CHLEAFXY;
		double[,] CHUCXY;
		double[,] CHV2XY;
		double[,] CHB2XY;
		double[,] ZLVLXY;
		
		public void UpdateMatrix()
		{
			for (int I = 0; I < NRows; I++) {
				for (int J = 0; J < NCols; J++) {
					GridCell cell = Cells[I, J];
					

					//#ifdef WRF_HYDRO
					//AD_CHANGE: Glacier cells can produce small negative subsurface runoff for mass balance.
					//       This will crash channel routing, so only pass along positive runoff.
					//soldrain[I, J] = Math.Max(cell.RUNSB * DT, 0);       //mm , underground runoff
					//INFXSRT[I, J] = RUNSF * DT;       //mm , surface runoff
//						//#}
//
//
//						// INPUT/OUTPUT
//
					TSK[I, J] = cell.TRAD;
					HFX[I, J] = cell.FSH;
					GRDFLX[I, J] = cell.SSOIL;
					SMSTAV[I, J] = 0.0;  // [maintained as Noah consistency]
					SMSTOT[I, J] = 0.0;  // [maintained as Noah consistency]
					
					//#ifdef WRF_HYDRO
//					ACCPRCP[I, J] += cell.PRCP * DT;
//					ACCECAN[I, J] += cell.ECAN * DT;
//					ACCETRAN[I, J] += cell.ETRAN * DT;
//					ACCEDIR[I, J] += cell.EDIR * DT;
					//#}
					if (cell.ALBEDO > -999) {
						ALBEDO[I, J] = cell.ALBEDO;
						if (XLAND[I, J] > 1)
							ALBEDO[I, J] = -1e36;
					}
					SNOWC[I, J] = cell.FSNO;
					for (int z = 1; z <= NSoil; z++) {
						SMOIS[I, z - 1, J] = cell.SMC[z];
						SH2O[I, z - 1, J] = cell.SH2O[z];
						TSLB[I, z - 1, J] = cell.STC[z];
						
					}
					SNOW[I, J] = cell.SNEQV;
					SNOWH[I, J] = cell.SNOWH;
					CANWAT[I, J] = cell.CANWAT;
					
					EMISS[I, J] = cell.EMISSI;
					QSFC[I, J] = cell.QSFC;

					ISNOWXY[I, J] = cell.ISNOW;
					TVXY[I, J] = cell.TV;
					TGXY[I, J] = cell.TG;
					CANLIQXY[I, J] = cell.CANLIQ;
					CANICEXY[I, J] = cell.CANICE;
					EAHXY[I, J] = cell.EAH;
					TAHXY[I, J] = cell.TAH;
					CMXY[I, J] = cell.CM;
					CHXY[I, J] = cell.CH;
					FWETXY[I, J] = cell.FWET;
					SNEQVOXY[I, J] = cell.SNEQVO;
					ALBOLDXY[I, J] = cell.ALBOLD;
					QSNOWXY[I, J] = cell.QSNOW;
					WSLAKEXY[I, J] = cell.WSLAKE;
					ZWTXY[I, J] = cell.ZWT;
					WAXY[I, J] = cell.WA;
					WTXY[I, J] = cell.WT;
					for (int z = 0; z < NSnow + NSoil; z++) {
						ZSNSOXY[I, z, J] = cell.ZSNSO.data[z];
					}
					for (int z = 0; z < NSnow; z++) {
						TSNOXY[I, z, J] = cell.STC.data[z];
						SNICEXY[I, z, J] = cell.SNICE.data[z];
						SNLIQXY[I, z, J] = cell.SNLIQ.data[z];
					}
					LFMASSXY[I, J] = cell.LFMASS;
					RTMASSXY[I, J] = cell.RTMASS;
					STMASSXY[I, J] = cell.STMASS;
					WOODXY[I, J] = cell.WOOD;
					STBLCPXY[I, J] = cell.STBLCP;
					FASTCPXY[I, J] = cell.FASTCP;
					LAIXY[I, J] = cell.LAI;
					XSAIXY[I, J] = cell.SAI;
					TAUSSXY[I, J] = cell.TAUSS;
					RAINBL[I, J] = cell.PRCP;
//					if(cell.PRCP!=0)
//						throw new Exception();

					// OUTPUT

					T2MVXY[I, J] = cell.T2MV;
					T2MBXY[I, J] = cell.T2MB;
					Q2MVXY[I, J] = cell.Q2V / (1.0 - cell.Q2V);  // specific humidity to mixing ratio
					Q2MBXY[I, J] = cell.Q2B / (1.0 - cell.Q2B); // consistent with registry def of Q2
					TRADXY[I, J] = cell.TRAD;
					NEEXY[I, J] = cell.NEE;
					GPPXY[I, J] = cell.GPP;
					NPPXY[I, J] = cell.NPP;
					FVEGXY[I, J] = cell.FVEGMP;
					VEGFRA[I, J] = cell.FVEG * 100; //这句并不保存
					RUNSFXY[I, J] = cell.RUNSRF;
					RUNSBXY[I, J] = cell.RUNSUB;
					ECANXY[I, J] = cell.ECAN;
					EDIRXY[I, J] = cell.EDIR;
					ETRANXY[I, J] = cell.ETRAN;
					FSAXY[I, J] = cell.FSA;
					FIRAXY[I, J] = cell.FIRA;
					APARXY[I, J] = cell.APAR;
					PSNXY[I, J] = cell.PSN;
					SAVXY[I, J] = cell.SAV;
					SAGXY[I, J] = cell.SAG;
					RSSUNXY[I, J] = cell.RSSUN;
					RSSHAXY[I, J] = cell.RSSHA;
					BGAPXY[I, J] = cell.BGAP;
					WGAPXY[I, J] = cell.WGAP;
					TGVXY[I, J] = cell.TGV;
					TGBXY[I, J] = cell.TGB;
					CHVXY[I, J] = cell.CHV;
					CHBXY[I, J] = cell.CHB;
					IRCXY[I, J] = cell.IRC;
					IRGXY[I, J] = cell.IRG;
					SHCXY[I, J] = cell.SHC;
					SHGXY[I, J] = cell.SHG;
					EVGXY[I, J] = cell.EVG;
					GHVXY[I, J] = cell.GHV;
					IRBXY[I, J] = cell.IRB;
					SHBXY[I, J] = cell.SHB;
					EVBXY[I, J] = cell.EVB;
					GHBXY[I, J] = cell.GHB;
					TRXY[I, J] = cell.TR;
					EVCXY[I, J] = cell.EVC;
					CHLEAFXY[I, J] = cell.CHLEAF;
					CHUCXY[I, J] = cell.CHUC;
					CHV2XY[I, J] = cell.CHV2;
					CHB2XY[I, J] = cell.CHB2;
					ZLVLXY[I, J] = cell.ZLVL; //2022年4月29日加，原Noah-MP没有此项。忽略该变量会导致restart后结果与正常执行有差距
					//RECHXY[I, J] = RECHXY[I, J] + cell.RECH * 1E3;//RECHARGE TO THE WATER TABLE
					//DEEPRECHXY[I, J] += cell.DEEPRECH;
					//SMCWTDXY[I, J] = cell.SMCWTD;
				}// ENDDO ILOOP                                                       // of I loop
			}//  ENDDO JLOOP
		}
		
		/// <summary>
		/// 输出Noah-MP状态
		/// </summary>
		public void WriteOutput()
		{
			
	
			if (outputPath == null)
				outputPath = forcePath + "/NoahMP_Output";
			if (!Directory.Exists(outputPath)) {
				Directory.CreateDirectory(outputPath);
			}
			string filename = outputPath + "/" + time0.Year + time0.Month.ToString("00") + time0.Day.ToString("00") + time0.Hour.ToString("00") + "_Output.nc";
			
			int ncid = 0;
			NC_LowAPI.nc_create(filename, 1, ref ncid);
			//int timeId=-1;
			int latId = -1;
			int lonId = -1;
			int SoilLayerId = -1;
			int SnowLayerId = -1;
			int SnowSoilLayerId = -1;
			//NC_LowAPI.nc_def_dim(ncid,"Time",NetCDF.NC_UNLIMITED,ref timeId);
			NC_LowAPI.nc_def_dim(ncid, "rows", NRows, ref latId);
			NC_LowAPI.nc_def_dim(ncid, "cols", NCols, ref lonId);
			NC_LowAPI.nc_def_dim(ncid, "SoilLayers", NSoil, ref SoilLayerId);
			NC_LowAPI.nc_def_dim(ncid, "SnowLayers", NSnow, ref SnowLayerId);
			NC_LowAPI.nc_def_dim(ncid, "SoilSnowLayers", NSnow + NSoil, ref SnowSoilLayerId);
			int[] SoilDims = new int[]{ SoilLayerId, latId, lonId };
			int[] SnowDims = new int[]{ SnowLayerId, latId, lonId };
			int[] SoilSnowDims = new int[]{ SnowSoilLayerId, latId, lonId };
			int[] Dims = new int[]{ latId, lonId };
			//分两步进行，outputNode为0时，实现变量的定义，为1时进行数据的填充
			for (outputMode = 0; outputMode < 2; outputMode++) {
				if (outputMode == 1) { //在进行第二步中数据的填充之前，先要关掉变量的定义阶段
					NC_LowAPI.nc_enddef(ncid);
				}
				
				add_to_output(ncid, FVEGXY, "FVEG", "Green Vegetation Fraction", "-");
				add_to_output(ncid, LAIXY, "LAI", "Leaf area index", "-");
				add_to_output(ncid, XSAIXY, "SAI", "Stem area index", "-");
				// Forcing
				add_to_output(ncid, SWDOWN, "SWFORC", "Shortwave forcing", "W m{-2}");
				//add_to_output(ncid,     , "COSZ"    , "Cosine of zenith angle"                    , "W m{-2}"               );
				add_to_output(ncid, GLW, "LWFORC", "Longwave forcing", "W m{-2}");
				add_to_output(ncid, RAINBL, "RAINRATE", "Precipitation rate", "kg m{-2} s{-1}");
				// Grid energy budget terms
				add_to_output(ncid, EMISS, "EMISS", "Grid emissivity", "");
				add_to_output(ncid, FSAXY, "FSA", "Total absorbed SW radiation", "W m{-2}");
				add_to_output(ncid, FIRAXY, "FIRA", "Total net LW radiation to atmosphere", "W m{-2}");
				add_to_output(ncid, GRDFLX, "GRDFLX", "Heat flux into the soil", "W m{-2}");
				add_to_output(ncid, HFX, "HFX", "Total sensible heat to atmosphere", "W m{-2}");
				add_to_output(ncid, LH, "LH", "Total latent heat to atmosphere", "W m{-2}");
				add_to_output(ncid, ECANXY, "ECAN", "Canopy water evaporation rate", "kg m{-2} s{-1}");
				add_to_output(ncid, ETRANXY, "ETRAN", "Transpiration rate", "kg m{-2} s{-1}");
				add_to_output(ncid, EDIRXY, "EDIR", "Direct from soil evaporation rate", "kg m{-2} s{-1}");
				add_to_output(ncid, ALBEDO, "ALBEDO", "Surface albedo", "-");
				// Grid water budget terms - in addition to above
				add_to_output(ncid, ACCUGDRUNOFF, "UGDRNOFF", "Accumulated underground runoff", "mm");
				add_to_output(ncid, ACCSFCRUNOFF, "SFCRNOFF", "Accumulatetd surface runoff", "mm");
				add_to_output(ncid, CANLIQXY, "CANLIQ", "Canopy liquid water content", "mm");
				add_to_output(ncid, CANICEXY, "CANICE", "Canopy ice water content", "mm");
				add_to_output(ncid, ZWTXY, "ZWT", "Depth to water table", "m");
				add_to_output(ncid, WAXY, "WA", "Water in aquifer", "kg m{-2}");
				add_to_output(ncid, WTXY, "WT", "Water in aquifer and saturated soil", "kg m{-2}");
				add_to_output(ncid, ACCPRCP, "ACCPRCP", "Accumulated precip", "mm");
				add_to_output(ncid, ACCECAN, "ACCECAN", "Accumulated canopy evap", "mm");
				add_to_output(ncid, ACCETRAN, "ACCETRAN", "Accumulated transpiration", "mm");
				add_to_output(ncid, ACCEDIR, "ACCEDIR", "Accumulated direct soil evap", "mm");
				add_to_output(ncid, RUNCOEF, "RUNCOEF", "Runoff coefficient (Acc. Runoff/Acc.Pcp)", "-");
				add_to_output(ncid, BIAS, "BIAS", "(Evapotranspiration + RunOff - Pcp)/Pcp", "-");

				// Additional needed to close the canopy energy budget
				add_to_output(ncid, SAVXY, "SAV", "Solar radiative heat flux absorbed by vegetation", "W m{-2}");
				add_to_output(ncid, TRXY, "TR", "Transpiration heat", "W m{-2}");
				add_to_output(ncid, EVCXY, "EVC", "Canopy evap heat", "W m{-2}");
				add_to_output(ncid, IRCXY, "IRC", "Canopy net LW rad", "W m{-2}");
				add_to_output(ncid, SHCXY, "SHC", "Canopy sensible heat", "W m{-2}");
				// Additional needed to close the under canopy ground energy budget
				add_to_output(ncid, IRGXY, "IRG", "Ground net LW rad", "W m{-2}");
				add_to_output(ncid, SHGXY, "SHG", "Ground sensible heat", "W m{-2}");
				add_to_output(ncid, EVGXY, "EVG", "Ground evap heat", "W m{-2}");
				add_to_output(ncid, GHVXY, "GHV", "Ground heat flux + to soil vegetated", "W m{-2}");
				// Needed to close the bare ground energy budget
				add_to_output(ncid, SAGXY, "SAG", "Solar radiative heat flux absorbed by ground", "W m{-2}");
				add_to_output(ncid, IRBXY, "IRB", "Net LW rad to atm bare", "W m{-2}");
				add_to_output(ncid, SHBXY, "SHB", "Sensible heat to atm bare", "W m{-2}");
				add_to_output(ncid, EVBXY, "EVB", "Evaporation heat to atm bare", "W m{-2}");
				add_to_output(ncid, GHBXY, "GHB", "Ground heat flux + to soil bare", "W m{-2}");
				// Above-soil temperatures
				add_to_output(ncid, TRADXY, "TRAD", "Surface radiative temperature", "K");
				add_to_output(ncid, TGXY, "TG", "Ground temperature", "K");
				add_to_output(ncid, TVXY, "TV", "Vegetation temperature", "K");
				add_to_output(ncid, TAHXY, "TAH", "Canopy air temperature", "K");
				add_to_output(ncid, TGVXY, "TGV", "Ground surface Temp vegetated", "K");
				add_to_output(ncid, TGBXY, "TGB", "Ground surface Temp bare", "K");
				add_to_output(ncid, T2MVXY, "T2MV", "2m Air Temp vegetated", "K");
				add_to_output(ncid, T2MBXY, "T2MB", "2m Air Temp bare", "K");
				add_to_output(ncid, Q2MVXY, "Q2MV", "2m mixing ratio vegetated", "kg/kg");
				add_to_output(ncid, Q2MBXY, "Q2MB", "2m mixing ratio bare", "kg/kg");
				add_to_output(ncid, EAHXY, "EAH", "Canopy air vapor pressure", "Pa");
				add_to_output(ncid, FWETXY, "FWET", "Wetted or snowed fraction of canopy", "fraction");
				
				// Snow and soil - 3D terms
				add_to_output_3d(ncid, ZSNSOXY, SoilSnowDims, NSnow + NSoil, "ZSNSO_SN", "Snow layer depths from snow surface", "m");
				add_to_output_3d(ncid, SNICEXY, SnowDims, NSnow, "SNICE", "Snow layer ice", "mm");
				add_to_output_3d(ncid, SNLIQXY, SnowDims, NSnow, "SNLIQ", "Snow layer liquid water", "mm");
				add_to_output_3d(ncid, TSLB, SoilDims, NSoil, "SOIL_T", "soil temperature", "K");
				add_to_output_3d(ncid, SMOIS, SoilDims, NSoil, "SOIL_M", "volumetric soil moisture", "m{3} m{-3}");
				add_to_output_3d(ncid, SH2O, SoilDims, NSoil, "SOIL_W", "liquid volumetric soil moisture", "m3 m-3");
				add_to_output_3d(ncid, TSNOXY, SnowDims, NSnow, "SNOW_T", "snow temperature", "K");
				// Snow - 2D terms
				add_to_output(ncid, SNOWH, "SNOWH", "Snow depth", "m");
				add_to_output(ncid, SNOW, "SNEQV", "Snow water equivalent", "kg m{-2}");
				add_to_output(ncid, QSNOWXY, "QSNOW", "Snowfall rate", "mm s{-1}");
				add_to_output(ncid, ISNOWXY, "ISNOW", "Number of snow layers", "count");
				add_to_output(ncid, SNOWC, "FSNO", "Snow-cover fraction on the ground", "");
				add_to_output(ncid, ACSNOW, "ACSNOW", "accumulated snow fall", "mm");
				add_to_output(ncid, ACSNOM, "ACSNOM", "accumulated melting water out of snow bottom", "mm");
				// Exchange coefficients
				add_to_output(ncid, CMXY, "CM", "Momentum drag coefficient", "");
				add_to_output(ncid, CHXY, "CH", "Sensible heat exchange coefficient", "");
				add_to_output(ncid, CHVXY, "CHV", "Exchange coefficient vegetated", "m s{-1}");
				add_to_output(ncid, CHBXY, "CHB", "Exchange coefficient bare", "m s{-1}");
				add_to_output(ncid, CHLEAFXY, "CHLEAF", "Exchange coefficient leaf", "m s{-1}");
				add_to_output(ncid, CHUCXY, "CHUC", "Exchange coefficient bare", "m s{-1}");
				add_to_output(ncid, CHV2XY, "CHV2", "Exchange coefficient 2-meter vegetated", "m s{-1}");
				add_to_output(ncid, CHB2XY, "CHB2", "Exchange coefficient 2-meter bare", "m s{-1}");
				// Carbon allocation model
				if (NoahMP.DVEG == 2 || NoahMP.DVEG == 5) {  //如果是动态植被启用了，则输出这些变量，否则不用输出
					add_to_output(ncid, LFMASSXY, "LFMASS", "Leaf mass", "g m{-2}");
					add_to_output(ncid, RTMASSXY, "RTMASS", "Mass of fine roots", "g m{-2}");
					add_to_output(ncid, STMASSXY, "STMASS", "Stem mass", "g m{-2}");
					add_to_output(ncid, WOODXY, "WOOD", "Mass of wood and woody roots", "g m{-2}");
					add_to_output(ncid, STBLCPXY, "STBLCP", "Stable carbon in deep soil", "g m{-2}");
					add_to_output(ncid, FASTCPXY, "FASTCP", "Short-lived carbon in shallow soil", "g m{-2}");
					add_to_output(ncid, NEEXY, "NEE", "Net ecosystem exchange", "g m{-2} s{-1} CO2");
					add_to_output(ncid, GPPXY, "GPP", "Net instantaneous assimilation", "g m{-2} s{-1} C");
					add_to_output(ncid, NPPXY, "NPP", "Net primary productivity", "g m{-2} s{-1} C");
					add_to_output(ncid, PSNXY, "PSN", "Total photosynthesis", "umol CO@ m{-2} s{-1}");
					add_to_output(ncid, APARXY, "APAR", "Photosynthesis active energy by canopy", "W m{-2}");
				}

				//        // Carbon allocation model
				//        if(RUNOFF_OPTION == 5); {
				//              add_to_output(ncid,SMCWTDXY   , "SMCWTD"   , "Leaf mass"                            , "g m{-2}"               );
				//              add_to_output(ncid,RECHXY     , "RECH"     , "Mass of fine roots"                   , "g m{-2}"               );
				//              add_to_output(ncid,QRFSXY     , "QRFS"     , "Stem mass"                            , "g m{-2}"               );
				//              add_to_output(ncid,QSPRINGSXY , "QSPRINGS" , "Mass of wood and woody roots"         , "g m{-2}"               );
				//              add_to_output(ncid,QSLATXY    , "QSLAT"    , "Stable carbon in deep soil"           , "g m{-2}"               );
				
			}
			NC_LowAPI.nc_close(ncid);
		}
		public void WriteRestart(string filename)
		{
			//string filename = time0.Year + time0.Month.ToString("00") + time0.Day.ToString("00") + time0.Hour.ToString("00") + "Output.nc";
			
			int ncid = 0;
			NC_LowAPI.nc_create(filename, 1, ref ncid);
			int TimeId = -9999;
			int DateLenId = -9999;
			int colsId = -9999;
			int rowsId = -9999;
			int dimid_dum = -9999;
			int dimid_layers = -9999;
			int dimid_snow_layers = -9999;
			int dimid_sosn_layers = -9999;
			NC_LowAPI.nc_def_dim(ncid, "Time", 0, ref TimeId); //0表示unlimited
			NC_LowAPI.nc_def_dim(ncid, "DateStrLen", 19, ref DateLenId);
			NC_LowAPI.nc_def_dim(ncid, "west_east", NCols, ref colsId);
			NC_LowAPI.nc_def_dim(ncid, "south_north", NRows, ref rowsId);
			NC_LowAPI.nc_def_dim(ncid, "west_east_stag", NCols + 1, ref dimid_dum);
			NC_LowAPI.nc_def_dim(ncid, "south_north_stag", NRows + 1, ref dimid_dum);
			NC_LowAPI.nc_def_dim(ncid, "soil_layers_stag", NSoil, ref dimid_layers);
			NC_LowAPI.nc_def_dim(ncid, "snow_layers", NSnow, ref dimid_snow_layers);
			NC_LowAPI.nc_def_dim(ncid, "sosn_layers", NSnow + NSoil, ref dimid_sosn_layers);
//			NC_LowAPI.nc_put_att_string(ncid, -1, "TITLE", 30, "RESTART FILE FROM Noah-MP-Opt");
			
			
			//分两步进行，outputNode为0时，实现变量的定义，为1时进行数据的填充
			for (outputMode = 0; outputMode < 2; outputMode++) {
				
				int[] SoilDims = new int[]{ rowsId, dimid_layers, colsId };
				int[] SnowDims = new int[]{ rowsId, dimid_snow_layers, colsId };
				int[] SoilSnowDims = new int[]{ rowsId, dimid_sosn_layers, colsId };
				int[] Dims = new int[]{ rowsId, colsId };
				
				add_to_restart_3d(ncid, TSLB, SoilDims, NSoil, "SOIL_T", "", "-");
				add_to_restart_3d(ncid, TSNOXY, SnowDims, NSnow, "SNOW_T", "", "-");
				add_to_restart_3d(ncid, SMOIS, SoilDims, NSoil, "SMC", "", "-");
				add_to_restart_3d(ncid, SH2O, SoilDims, NSoil, "SH2O", "", "-");
				add_to_restart_3d(ncid, ZSNSOXY, SoilSnowDims, NSnow + NSoil, "ZSNSO", "", "-");
				add_to_restart_3d(ncid, SNICEXY, SnowDims, NSnow, "SNICE", "", "-");
				add_to_restart_3d(ncid, SNLIQXY, SnowDims, NSnow, "SNLIQ", "", "");
				add_to_restart(ncid, QSNOWXY, Dims, "QSNOW", "", "");
				add_to_restart(ncid, FWETXY, Dims, "FWET", "", "");
				add_to_restart(ncid, SNEQVOXY, Dims, "SNEQVO", "", "");
				add_to_restart(ncid, EAHXY, Dims, "EAH", "", "");
				add_to_restart(ncid, TAHXY, Dims, "TAH", "", "");
				add_to_restart(ncid, ALBOLDXY, Dims, "ALBOLD", "", "");
				add_to_restart(ncid, CMXY, Dims, "CM", "", "");
				add_to_restart(ncid, CHXY, Dims, "CH", "", "");
				add_to_restart(ncid, ISNOWXY, Dims, "ISNOW", "", "");
				add_to_restart(ncid, CANLIQXY, Dims, "CANLIQ", "", "");
				add_to_restart(ncid, CANICEXY, Dims, "CANICE", "", "");
				add_to_restart(ncid, SNOW, Dims, "SNEQV", "", "");
				add_to_restart(ncid, SNOWH, Dims, "SNOWH", "", "");
				add_to_restart(ncid, TVXY, Dims, "TV", "", "");
				add_to_restart(ncid, TGXY, Dims, "TG", "", "");
				add_to_restart(ncid, ZWTXY, Dims, "ZWT", "", "");
				add_to_restart(ncid, WAXY, Dims, "WA", "", "");
				add_to_restart(ncid, WTXY, Dims, "WT", "", "");
				add_to_restart(ncid, WSLAKEXY, Dims, "WSLAKE", "", "");
				add_to_restart(ncid, LFMASSXY, Dims, "LFMASS", "", "");
				add_to_restart(ncid, RTMASSXY, Dims, "RTMASS", "", "");
				add_to_restart(ncid, STMASSXY, Dims, "STMASS", "", "");
				add_to_restart(ncid, WOODXY, Dims, "WOOD", "", "");
				add_to_restart(ncid, STBLCPXY, Dims, "STBLCP", "", "");
				add_to_restart(ncid, FASTCPXY, Dims, "FASTCP", "", "");
				add_to_restart(ncid, LAIXY, Dims, "LAI", "", "");
				add_to_restart(ncid, XSAIXY, Dims, "SAI", "", "");
				add_to_restart(ncid, VEGFRA, Dims, "FPAR", "", "");
				add_to_restart(ncid, GVFMIN, Dims, "GVFMIN", "", "");
				add_to_restart(ncid, GVFMAX, Dims, "GVFMAX", "", "");
				add_to_restart(ncid, SHDMAX, Dims, "SHDMAX", "", "");
				add_to_restart(ncid, ACSNOM, Dims, "ACMELT", "", "");
				add_to_restart(ncid, ACSNOW, Dims, "ACSNOW", "", "");
				add_to_restart(ncid, TAUSSXY, Dims, "TAUSS", "", "");
				add_to_restart(ncid, QSFC, Dims, "QSFC", "", "");
				add_to_restart(ncid, ACCSFCRUNOFF, Dims, "SFCRUNOFF", "", "");
				add_to_restart(ncid, ACCUGDRUNOFF, Dims, "UGDRUNOFF", "", "");
				add_to_restart(ncid, ACCPRCP, Dims, "ACCPRCP", "", "");
				add_to_restart(ncid, ACCECAN, Dims, "ACCECAN", "", "");
				add_to_restart(ncid, ACCEDIR, Dims, "ACCEDIR", "", "");
				add_to_restart(ncid, ACCETRAN, Dims, "ACCETRAN", "", "");
				if (NoahMP.OPT_RUN == 5) {
					add_to_restart_3d(ncid, SMOISEQ, SoilDims, NSoil, "SMOISEQ", "", "");
					add_to_restart(ncid, AREAXY, Dims, "AREAXY", "", "");
					add_to_restart(ncid, SMCWTDXY, Dims, "SMCWTDXY", "", "");
					add_to_restart(ncid, DEEPRECHXY, Dims, "DEEPRECHXY", "", "");
					add_to_restart(ncid, QSLATXY, Dims, "QSLATXY", "", "");
					add_to_restart(ncid, QRFSXY, Dims, "QRFSXY", "", "");
					add_to_restart(ncid, QSPRINGSXY, Dims, "QSPRINGSXY", "", "");
					add_to_restart(ncid, RECHXY, Dims, "RECHXY", "", "");
					add_to_restart(ncid, QRFXY, Dims, "QRFXY", "", "");
					add_to_restart(ncid, QSPRINGXY, Dims, "QSPRINGXY", "", "");
					add_to_restart(ncid, FDEPTHXY, Dims, "FDEPTHXY", "", "");
					add_to_restart(ncid, RIVERCONDXY, Dims, "RIVERCONDXY", "", "");
					add_to_restart(ncid, RIVERBEDXY, Dims, "RIVERBEDXY", "", "");
					add_to_restart(ncid, EQZWT, Dims, "EQZWT", "", "");
					add_to_restart(ncid, PEXPXY, Dims, "PEXPXY", "", "");
				}
				add_to_restart(ncid, ZLVLXY, Dims, "ZLVLXY", "", "");//2022年4月29日加，原Noah-MP没有此项。忽略该变量会导致restart后结果与正常执行有差距
				if (outputMode == 0) { //在进行第二步中数据的填充之前，先要关掉变量的定义阶段
					NC_LowAPI.nc_enddef(ncid);
				} else {
					NC_LowAPI.nc_close(ncid);
				}
			}
			//记录下上一次写入盘的文件名
			StreamWriter sw = new StreamWriter("LastRestart.txt");
			sw.WriteLine(filename);
			sw.Close();
		}
	}
}


