/*
 * Created by SharpDevelop.
 * User: Administrator
 * Date: 2021/3/14
 * Time: 20:40
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;
using System.Collections.Generic;
using System.IO;
using System.Threading.Tasks;

namespace NoahMP
{
	public struct RouteFlux
	{
		public double influx;
		public double outflux;
		/// <summary>
		/// 容量
		/// </summary>
		public double capacity;
		/// <summary>
		/// 记录接收上游来水的个数
		/// </summary>
		public int nUpperStream;
		public bool mark;
	}
	public class StreamInfo
	{
		/// <summary>
		/// 在河网中的id号
		/// </summary>
		public long hydroId;
		/// <summary>
		/// 在河网中下一段的id号
		/// </summary>
		public long nextId;
		/// <summary>
		/// 河段的汇流区面积
		/// </summary>
		public double area;
		/// <summary>
		/// 行号，在runoff网格中的行号(指粗格点)
		/// </summary>
		public int row;
		/// <summary>
		/// 列号,在runoff网格中的列号(指粗格点)
		/// </summary>
		public int col;
		/// <summary>
		/// 下一河段在数组中的序号，与前面的id号不同。数组序号用于汇流演算时的快速访问
		/// </summary>
		public int nextindex;
		/// <summary>
		/// 上游入流的个数
		/// </summary>
		public int upperCount;
		public double length;
		/// <summary>
		/// 流量
		/// </summary>
		public double outQ;
		public double capacity;
		/// <summary>
		/// 跨计算核的关键节点，在进入边界处，需要特殊处理
		/// </summary>
		public bool IsConjuction;
	}
	/// <summary>
	/// Description of RiverRouting.
	/// </summary>
	public class RiverRouting
	{
		public static int threads = 4;
		/// <summary>
		/// wave celerity, unit:km/h
		/// </summary>
		public double c0 = 1;
		
		/// <summary>
		/// 0.35
		/// </summary>
		public double lambd = 0.35;
		//
		/// <summary>
		/// 900s  //The calculation time step,15 min is default. unit: s
		/// </summary>
		public double TimeStep = 0;
		/// <summary>
		///  10800s
		/// </summary>
		public double DataStep = 0;
		//
		public double X = 0.4;
		//0.3
		public double AreaFactor = 0;
		RouteFlux[] preQ;
		RouteFlux[] curQ;
		StreamInfo[] ChanInfo;
		StreamWriter fp = null;
		string path;
		//		public double[,] runoff0;
		//		public double[,] runoff1;
		StreamInfo[] order;
		SortedList<int,int> slist = new SortedList<int, int>();
		/// <summary>
		/// 收集每块中的交叉河段
		/// </summary>
		List<int>[] slistIndex;
		public RiverRouting()
		{
			//			int count,ti,hi,i;
			//global preQ,curQ,StrInfo,fp
			order = ReadChannelOrder("ChannelOrder.txt");
			ReadParams();
			int count = order.Length;
			preQ = new RouteFlux[count];
			curQ = new RouteFlux[count];
			ChanInfo = new StreamInfo[count];
			//preQ= <RouteFlux*> malloc(count*sizeof(RouteFlux))
			//curQ= <RouteFlux*> malloc(count*sizeof(RouteFlux))
			//StrInfo= <StreamInfo*> malloc(count*sizeof(StreamInfo))    	
			InitOrderlist(order);	
			slist.Clear();
			//threads = 16;
			slistIndex = new List<int>[threads];
			int num = count / threads + 1;
			
			for (int k = 0; k < threads; k++) {
				slistIndex[k] = new List<int>();
				int start = k * num;
				for (int i = start; i < start + num; i++) {
					if (i >= ChanInfo.Length)
						break;
					int nextindex = ChanInfo[i].nextindex;
					//如果当前河段的下游河段位于当前块之外，则将下游河段标记为跨块河段，将跨块河段收集到列表中
					if (nextindex >= start + num && nextindex != -1) {						
						
						ChanInfo[nextindex].IsConjuction = true;	
						if (!slist.ContainsKey(nextindex)) {
							slistIndex[k].Add(nextindex);
						}
						slist[nextindex] = 1;
					}
					//如果列表中已包含了当前河段(即当前河段是跨块河段)，其下游河段也要被收入列表中
					if (slist.ContainsKey(i) && nextindex != -1) {
						
						if (!slist.ContainsKey(nextindex)) {
							slistIndex[k].Add(nextindex);
						}
						slist[nextindex] = 1;
					}
				}
			}
			Console.WriteLine(threads + " " + num + " " + slist.Count + " " + (num + slist.Count));
					
		}
		
		public StreamInfo[] ReadChannelOrder(string filename)
		{
			List<StreamInfo> mList = new List<StreamInfo>();
			StreamReader sr = new StreamReader(filename);
			sr.ReadLine();
			//读入河网
			while (!sr.EndOfStream) {
				string line = sr.ReadLine();	
				StreamInfo info = new StreamInfo();	
				info.hydroId = -1;
				if (line.Contains("exceeds")) {
					mList.Add(info);
					continue;
				}
				
				string[] strs = line.Split(new char[]{ ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
				//StreamInfo info = new StreamInfo();
				info.hydroId = (long)Convert.ToDouble(strs[0]);
				info.nextId = (long)Convert.ToDouble(strs[1]);
				info.area = Convert.ToDouble(strs[2]);
				info.row = Convert.ToInt32(strs[3]);
				info.col = Convert.ToInt32(strs[4]);
				info.length = Convert.ToDouble(strs[5]);
				mList.Add(info);
			}
			sr.Close();
			//准备河段在数组中的序号，用于模拟过程中快速访问
			for (int i = 0; i < mList.Count; i++) {
				StreamInfo info = mList[i];
				if (info == null) {
					continue;
				}
				int index = -1;
				for (int j = 0; j < mList.Count; j++) {
					StreamInfo info2 = mList[j];
					if (info2 != null && info2.hydroId == info.nextId) {
						index = j;
						info2.upperCount += 1;
						break;
					}
				}
				info.nextindex = index;
			}
			
			return mList.ToArray();
		}
		public void ReadParams()
		{
			StreamReader sr = new StreamReader("namelist");
			while (!sr.EndOfStream) {
				string line = sr.ReadLine();    
				line = line.Split(new char[]{ '!' }, StringSplitOptions.RemoveEmptyEntries)[0];
				string[] mline = line.Split(new char[]{ '=' }, StringSplitOptions.RemoveEmptyEntries);
				// print line
				if (mline[0] == "waveCelerity") {
					c0 = Convert.ToDouble(mline[1]);
					continue;
				}
				if (mline[0] == "lambda") {
					lambd = Convert.ToDouble(mline[1]);
					continue;
				}
				if (mline[0] == "X") {
					X = Convert.ToDouble(mline[1]);
					continue;
				}
				if (mline[0] == "DataStep") {
					DataStep = Convert.ToDouble(mline[1]);
					continue;
				}
				if (mline[0] == "RoutingStep") {
					TimeStep = Convert.ToDouble(mline[1]);
					continue;
				}
				if (mline[0] == "AreaPerCell") {
					AreaFactor = Convert.ToDouble(mline[1]);
					continue;
				}
			}
		}
		public void InitOrderlist(StreamInfo[] orderlist)
		{
			ChanInfo = orderlist;
			for (int i = 0; i < orderlist.Length; i++) {
				preQ[i].influx = 0;
				preQ[i].outflux = 0;
				preQ[i].capacity = 0;
				curQ[i].influx = 0;
				curQ[i].outflux = 0;
				curQ[i].capacity = 0;
			}
		}
		public double Muskingum(double Q11, double Q12, double Q21, double dGridQ, double x, double k, double dT)
		{
			double denominator, C1, C2, C3, C4, Q22;
			denominator = k * (1 - x) + dT / 2;
			C1 = (dT / 2 + x * k) / denominator;
			C2 = (dT / 2 - x * k) / denominator;
			C3 = 1 - C1 - C2;  //(k * (1 - x) - dT/2 ) / denominator;
			C4 = dT / denominator;   //dGridQ * dT * Lm / denominator;
			//the first subscript represent position, the second subscript represent time 
			//计算Cunge法的有限差分值: 流量=C1*上段面t时流量+C2*上段面t+1时流量+C3*下段面t时流量
			//当某段河流上一时刻入流和出流均为0且当前时刻有入流时，且C2<0时，则会计算出现负值，
			//这种情况先将当前入流的水分配一半至出流？
			Q22 = C1 * Q11 + C2 * Q12 + C3 * Q21 + C4 * dGridQ;
			return Q22;
		}
		public void Routing(double[,] runoff, int hi)
		{
			//area,cc,length;
			RouteFlux[,] temp;     	
			//global preQ,curQ,StrInfo,fp
			int mcount = (int)(DataStep / TimeStep);
			
			for (int i = 0; i < ChanInfo.Length; i++) {
				if (ChanInfo[i] == null)
					continue;
				int row = ChanInfo[i].row;
				int col = ChanInfo[i].col;
				double area = ChanInfo[i].area * AreaFactor;
				double length = ChanInfo[i].length;
				//c0 = ChanInfo[i].celerity(preQ[i].outflux);
				double cc = c0 * 1000.0 / 3600;   // change unit to m/s
				double K = lambd * length / cc;  //lambd=0.35, as estimated by Cedric David 2011
				//double K = length / cc;
				//print X,K,TimeStep,DataStep,c0,lambd
				int nextIndex = ChanInfo[i].nextindex;
				//if runoff[row,col]<0:
				//    runoff[row,col]=0
				double tempQ = runoff[row, col] * area * 1000 / 12;   //*(TimeStep/DataStep)   //unit: m^3 per calculation step

				//计算Cunge法的有限差分值: 流量=C1*时段前上段面流量+C2*时段后上段面流量+C3*时段前下段面流量
				double outQ = Muskingum(preQ[i].influx, curQ[i].influx, preQ[i].outflux, tempQ, X, K, TimeStep);
				if (double.IsNaN(outQ))
					throw new Exception();
				curQ[i].outflux = outQ;
				ChanInfo[i].outQ = outQ;
				ChanInfo[i].capacity -= outQ;
				if (nextIndex != -1) {
					curQ[nextIndex].influx += outQ;
					ChanInfo[nextIndex].capacity += outQ;
				}
				
			}
		}
		
		public int ProcessChn(double[,] runoff, int i)
		{
			if (ChanInfo[i].hydroId < 0)
				return -1;
			int row = ChanInfo[i].row;
			int col = ChanInfo[i].col;
			double area = ChanInfo[i].area * AreaFactor;
			double length = ChanInfo[i].length;
			double cc = c0 * 1000.0 / 3600;   // change unit to m/s
			double K = lambd * length / cc;
			//print X,K,TimeStep,DataStep,c0,lambd
			int nextIndex = ChanInfo[i].nextindex;
			double tempQ = runoff[row, col] * area * 1000 / 12;   //*(TimeStep/DataStep)   //unit: m^3 per calculation step

			//计算Cunge法的有限差分值: 流量=C1*时段前上段面流量+C2*时段后上段面流量+C3*时段前下段面流量
			double outQ = Muskingum(preQ[i].influx, curQ[i].influx, preQ[i].outflux, tempQ, X, K, TimeStep);
			curQ[i].outflux = outQ;
			curQ[i].capacity -= outQ;
			if (nextIndex != -1) {
				curQ[nextIndex].influx += outQ;
				curQ[nextIndex].capacity += outQ;
				curQ[nextIndex].nUpperStream += 1;
			}
			curQ[i].mark = true;
			ChanInfo[i].outQ = outQ;
			return nextIndex;
		}
		public int ProcessChn(double[,] runoff, int i, bool flag)
		{
			if (ChanInfo[i].hydroId < 0)
				return -1;
			int row = ChanInfo[i].row;
			int col = ChanInfo[i].col;
			double area = ChanInfo[i].area * AreaFactor;
			double length = ChanInfo[i].length;
			
			//c0 = ChanInfo[i].celerity(preQ[i].outflux);
			double cc = c0 * 1000.0 / 3600;   // change unit to m/s
			double K = lambd * length / cc;  //lambd=0.35, as estimated by Cedric David 2011
			//double K = length / cc;
			
			//print X,K,TimeStep,DataStep,c0,lambd
			int nextIndex = ChanInfo[i].nextindex;
			double tempQ = runoff[row, col] * area * 1000 / 12;   //*(TimeStep/DataStep)   //unit: m^3 per calculation step

			//计算Cunge法的有限差分值: 流量=C1*时段前上段面流量+C2*时段后上段面流量+C3*时段前下段面流量
			double outQ = 0;
			if (flag) {
				outQ = Muskingum(preQ[i].influx, preQ[i].influx, preQ[i].outflux, tempQ, X, K, TimeStep);
			} else {
				outQ = Muskingum(preQ[i].influx, curQ[i].influx, preQ[i].outflux, tempQ, X, K, TimeStep);
			}
			curQ[i].outflux = outQ;
			ChanInfo[i].capacity -= outQ;
			if (nextIndex != -1) {
				curQ[nextIndex].influx += outQ;
				ChanInfo[nextIndex].capacity += outQ;
				curQ[nextIndex].nUpperStream += 1;
			}
			curQ[i].mark = true;
			ChanInfo[i].outQ = outQ;
			return nextIndex;
		}
		/// <summary>
		/// 用于修正跨界河段
		/// </summary>
		/// <param name="runoff"></param>
		/// <param name="i"></param>
		/// <param name="i">下一块的河段起始编号</param>
		/// <returns></returns>
		public int ProcessChnStep2(double[,] runoff, int i,int nextStart)
		{
			if (ChanInfo[i].hydroId < 0)
				return -1;
			int row = ChanInfo[i].row;
			int col = ChanInfo[i].col;
			double area = ChanInfo[i].area * AreaFactor;
			double length = ChanInfo[i].length;
			//c0 = ChanInfo[i].celerity(preQ[i].outflux);
			double cc = c0 * 1000.0 / 3600;   // change unit to m/s
			double K = lambd * length / cc;  //lambd=0.35, as estimated by Cedric David 2011
			//double K = length / cc;
			
			//print X,K,TimeStep,DataStep,c0,lambd
			int nextIndex = ChanInfo[i].nextindex;
			double tempQ = runoff[row, col] * area * 1000 / 12;   //*(TimeStep/DataStep)   //unit: m^3 per calculation step

			//计算Cunge法的有限差分值: 流量=C1*时段前上段面流量+C2*时段后上段面流量+C3*时段前下段面流量
			double outQ = Muskingum(preQ[i].influx, curQ[i].influx, preQ[i].outflux, tempQ, X, K, TimeStep);
			double delta = outQ - curQ[i].outflux;
//			Console.WriteLine(outQ.ToString("0.0000")+" "+delta.ToString("0.0000")+" "+(delta/outQ).ToString("0.0000"));
//			if(delta>1e-10)
//			Thread.Sleep(3000);
			curQ[i].outflux = outQ;
			
			if (nextIndex != -1 && nextIndex<nextStart) { 
				curQ[nextIndex].influx += delta;
			}
			//curQ[i].mark = true;
			ChanInfo[i].outQ = outQ;
			return nextIndex;
		}
		public void RoutingParallel(double[,] runoff, int hi)
		{
			RouteFlux[,] temp;     				
			int nThread = NoahMP.threads;
			int count = ChanInfo.Length / nThread + 1;
			Parallel.For(0, nThread, k => {
				//for (int k = 0; k < nThread; k++) {
				int start = k * count;
				//StreamWriter sw=new StreamWriter("channel_marks.txt"+k);
				//在本块内部时要按序列的存储顺序处理
				for (int i = start; i < start + count; i++) {
					if (i >= ChanInfo.Length)
						break;
					if (slist.ContainsKey(i)) {
						//sw.WriteLine(ChanInfo[i].hydroId+" -1");
						continue;
					}
					//sw.WriteLine(ChanInfo[i].hydroId + " " + k);
					ProcessChn(runoff, i);					
				}
				//sw.Close();
			});
			
			for (int ind = 0; ind < slist.Count; ind++) {
				int i = slist.Keys[ind];
				ProcessChn(runoff, i);
			}
			
		}
		/// <summary>
		/// RAPID法
		/// </summary>
		/// <param name="runoff"></param>
		/// <param name="hi"></param>
		public void RoutingParallel3(double[,] runoff, int hi)
		{
			RouteFlux[,] temp;     				
			int nThread = threads;
			int count = ChanInfo.Length / nThread + 1;
			ParallelOptions option = new ParallelOptions();
			option.MaxDegreeOfParallelism = 4;
			Parallel.For(0, nThread, option, k => {
				int start = k * count;				
				//在本块内部时要按序列的存储顺序处理
				for (int i = start; i < start + count; i++) {
					if (i >= ChanInfo.Length) {
						break;
					}					
					if (ChanInfo[i].IsConjuction) {
						ProcessChn(runoff, i, true);
					} else {
						ProcessChn(runoff, i);	
					}					
				}
			});
			
			Parallel.For(0, nThread, option, k => {
			//for (int k = 0; k < nThread; k++) {
			    int start = k * count;
			    int nextStart=start+count;
				List<int> setn = slistIndex[k];
				for (int index = 0; index < setn.Count; index++) {			        	
					int i = setn[index];
					ProcessChnStep2(runoff, i,nextStart);
					//ProcessChn(runoff, i);
				}
			//}
			});
		}
		
		/// <summary>
		/// 
		/// </summary>
		/// <param name="path"></param>
		/// <param name="year"></param>
		/// <param name="month"></param>
		/// <param name="day"></param>
		/// <param name="hour"></param>
		/// <param name="ntimes"></param>
		/// <param name="SFCRNOFF"></param>
		/// <param name="UGDRNOFF"></param>
		public void RunRouting(DateTime time0, double[,] TOTALRUNOFF)
		{
			
//			if (month == 1 && day == 1 && hour == 0) {
//				if(fp!=null){
//					fp.Close();
//				}
//				//print year
//				fp = new StreamWriter("channel_" + year + ".txt");
//			}			
			
			//每15分钟一次，三小时共12次
			int count = order.Length;
			for (int hi = 0; hi < 12; hi++) {
				//Routing(TOTALRUNOFF, hi);
				RoutingParallel3(TOTALRUNOFF, hi);
				RouteFlux[] temp = preQ;
				preQ = curQ;
				curQ = temp;
				for (int i = 0; i < count; i++) {
//					if (ChanInfo[i].hydroId != -1)// && !preQ[i].mark)
//						throw new Exception();
					curQ[i].influx = 0;
					curQ[i].outflux = 0;
					curQ[i].nUpperStream = 0;
					curQ[i].mark = false;
				}
			}
			//保存
			int year = time0.Year;
			int month = time0.Month;
			int day = time0.Day;
			int hour = time0.Hour;
			if ((month == 1 && day == 1 && hour == 0) || fp == null) {
				if (Driver.outputPath == null)
					Driver.outputPath = Driver.forcePath + "/NoahMP_Output";
				if (!Directory.Exists(Driver.outputPath)) {
					Directory.CreateDirectory(Driver.outputPath);
				}
				string fileName = Driver.outputPath + "/channel_" + year.ToString("0000") + month.ToString("00") + day.ToString("00") + hour.ToString("00") + ".txt";
				if(fp!=null)
					fp.Close();
				fp = new StreamWriter(fileName);
			}
			fp.Write(year.ToString("0000") + month.ToString("00") + day.ToString("00") + hour.ToString("00") + "\t");
			for (int i = 0; i < count; i++) {
				fp.Write((ChanInfo[i].outQ / TimeStep).ToString("0.00") + "\t");
			}
			fp.WriteLine();
			//fp.Close();
				
		}
		public void WriteStatus(string fileName)
		{
			FileStream stream = new FileStream(fileName, FileMode.Create);
			BinaryWriter sw = new BinaryWriter(stream);
			sw.Write(preQ.Length);
			for (int i = 0; i < preQ.Length; i++) {
				sw.Write(preQ[i].influx);
				sw.Write(preQ[i].outflux);
				sw.Write(preQ[i].capacity);
			}
			sw.Close();
			stream.Close();
		}
		public void ReadStatus(string fileName)
		{
			FileStream stream = new FileStream(fileName, FileMode.Open);
			BinaryReader sw = new BinaryReader(stream);
			int count = sw.ReadInt32();
			for (int i = 0; i < count; i++) {
				preQ[i].influx = sw.ReadDouble();//preQ[i].influx);
				preQ[i].outflux = sw.ReadDouble();//preQ[i].outflux);
				preQ[i].capacity = sw.ReadDouble();//(preQ[i].capacity);
			}
			sw.Close();
			stream.Close();
		}
	}

	
}
