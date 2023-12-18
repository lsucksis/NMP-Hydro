/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2018/1/14
 * Time: 19:16
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;
using System.Threading;

namespace NoahMP
{
	/// <summary>
	/// Description of NoahMP_Cell.
	/// </summary>
	public class GridCellReal:GridCell
	{	
		
		public override void SurfaceFlux()
		{
			double COSZ = Driver.CALC_DECLIN(Driver.time0, this.LATITUDE, this.LONGITUDE);
			int YEARLEN = getYearLen(Driver.time0.Year);
			if(TG<100)
				throw new Exception("");
			if(this.SHDMAX>1 || this.SHDMAX<0.01 || this.SHDMIN>1)
				throw new Exception("");
		
			NoahMP3.NoahMP_SFLX(this, -1, -1, this.LATITUDE, YEARLEN, -1, COSZ,
				Driver.DT, Driver.DX, DZ8W, Driver.NSoil, ZSOIL, Driver.NSnow,
				FVEG, SHDMAX, VEGTYP, Driver.vegparams.ISURBAN, ICE, IST,
				ISC, SMCEQ,
				IZ0TLND, SFCTMP, SFCPRS, PSFC, UU, VV, Q2,
				QC, SOLDN, LWDN, PRCP, TBOT, CO2AIR,
				O2AIR, FOLN, FICEOLD, PBLH, ZLVL,
				ref ALBOLD,ref SNEQVO,
				STC, SH2O, SMC,ref TAH,ref EAH,ref FWET,
				ref CANLIQ, ref CANICE, ref TV, ref TG, ref QSFC,ref QSNOW,
				ref ISNOW, ZSNSO, ref SNOWH,ref SNEQV, SNICE, SNLIQ,
				ref ZWT, ref WA, ref WT, ref WSLAKE, ref LFMASS, ref RTMASS,
				ref STMASS, ref WOOD, ref STBLCP, ref FASTCP, ref LAI, ref SAI,
				ref CM, ref CH, ref TAUSS,
				ref SMCWTD, ref DEEPRECH, ref RECH,
				out FSA, out FSR, out FIRA, out FSH, out SSOIL, out FCEV,
				out FGEV, out FCTR, out ECAN, out ETRAN, out EDIR, out TRAD,
				out TGB, out TGV, out T2MV, out T2MB, out Q2V, out Q2B,
				out RUNSRF, out RUNSUB, out APAR, out PSN, out SAV, out SAG,
				out FSNO, out NEE, out GPP, out NPP, out FVEGMP, out ALBEDO,
				out QSNBOT, out PONDING, out PONDING1, out PONDING2, out RSSUN, out RSSHA,
				out BGAP, out WGAP, out CHV, out CHB, out EMISSI,
				out SHG, out SHC, out SHB, out EVG, out EVB, out GHV,
				out GHB, out IRG, out IRC, out IRB, out TR, out EVC,
				out  CHLEAF, out CHUC, out CHV2, out CHB2, out FPICE,
				ref  SFCHEADRT);
			
		}	
		
	}
}
